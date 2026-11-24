#!/usr/bin/env python3

import os
import re
import json
import argparse
import requests
import pandas as pd

DATABASES = [
    "GTEx_Tissues_V8_2023",
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "MSigDB_Hallmark_2020",
]

ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Functional enrichment analysis of gene modules and tumors"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV containing module, gene and disease information"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    return parser.parse_args()


def read_input(path):
    df = pd.read_csv(path)

    required_columns = {
        "module",
        "gene",
        "disease",
        "mean_expression",
        "pct_expressing",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {', '.join(missing)}"
        )
    df = df.dropna(subset=["gene", "disease"])

    df["gene"] = df["gene"].astype(str).str.strip()
    df["disease"] = df["disease"].astype(str).str.strip()

    df["mean_expression"] = pd.to_numeric(
        df["mean_expression"],
        errors="coerce"
    )

    df["pct_expressing"] = pd.to_numeric(
        df["pct_expressing"],
        errors="coerce"
    )

    return df


def read_modules(df):
    module_df = df.dropna(subset=["module"]).copy()

    modules = (
        module_df.groupby("module")["gene"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    return modules


def read_tumors(df):
    # Keep only genes with evidence of expression
    expressed = df[
        (df["mean_expression"] > 0) |
        (df["pct_expressing"] > 0)
    ].copy()

    tumors = (
        expressed.groupby("disease")["gene"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    return tumors


def upload_gene_list(genes, description):
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, description)
    }

    r = requests.post(
        ADD_URL,
        files=payload,
        timeout=120
    )
    r.raise_for_status()

    return json.loads(r.text)["userListId"]


def enrich(user_list_id, database, n_input_genes):
    url = (
        f"{ENRICH_URL}?userListId={user_list_id}"
        f"&backgroundType={database}"
    )

    r = requests.get(url, timeout=120)
    r.raise_for_status()

    data = json.loads(r.text).get(database, [])

    rows = []

    for x in data:
        matched_genes = x[5]
        rows.append({
            "rank": x[0],
            "term": x[1],
            "p_value": x[2],
            "odds_ratio": x[3],
            "combined_score": x[4],
            "genes": ";".join(matched_genes),
            "n_input_genes": n_input_genes,
            "n_matched_genes": len(matched_genes),
            "adjusted_p_value": x[6],
        })

    return pd.DataFrame(rows)


def run_enrichment(groups, output_dir, group_type):

    os.makedirs(output_dir, exist_ok=True)

    summary = []

    for group, genes in groups.items():

        if len(genes) == 0:
            continue
        safe_name = re.sub(
            r"[^a-zA-Z0-9_-]+",
            "_",
            str(group).strip()
        ).strip("_").lower()

        group_dir = os.path.join(
            output_dir,
            safe_name
        )

        os.makedirs(group_dir, exist_ok=True)

        user_list = upload_gene_list(
            genes,
            f"{group_type}: {group}"
        )

        for db in DATABASES:
            df = enrich(
                user_list,
                db,
                len(genes)
            )
            df.to_csv(
                os.path.join(group_dir, f"{db}.csv"),
                index=False
            )

            if not df.empty:
                top = df.iloc[0]
                summary.append({
                    group_type: group,
                    "database": db,
                    "term": top["term"],
                    "adjusted_p_value": top["adjusted_p_value"],
                    "combined_score": top["combined_score"],
                    "n_input_genes": top["n_input_genes"],
                    "n_matched_genes": top["n_matched_genes"],
                    "genes": top["genes"],
                })

    return pd.DataFrame(summary)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    df = read_input(args.input)

    # ==========================================
    # ENRICHMENT BY MODULE
    # ==========================================

    modules = read_modules(df)
    module_summary = run_enrichment(
        modules,
        os.path.join(args.outdir, "by_module"),
        "module"
    )
    module_summary.to_csv(
        os.path.join(args.outdir, "summary.csv"),
        index=False
    )

    # ==========================================
    # ENRICHMENT BY TUMOR
    # ==========================================

    tumors = read_tumors(df)
    tumor_summary = run_enrichment(
        tumors,
        os.path.join(args.outdir, "by_tumor"),
        "tumor"
    )
    tumor_summary.to_csv(
        os.path.join(args.outdir, "summary_by_tumor.csv"),
        index=False
    )

if __name__ == "__main__":
    main()
