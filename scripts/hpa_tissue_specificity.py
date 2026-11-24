
#!/usr/bin/env python3

import argparse
import json
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
from pathlib import Path

import pandas as pd
import requests


# CONFIGURATION
HPA_API_URL = "https://www.proteinatlas.org/api/search_download.php"

HPA_COLUMNS = (
    "g,eg,rnats,rnatd,rnatss,rnatsm,"
    "rnascs,rnascd,rnascss,rnascsm"
)

OUTPUT_COLUMNS = [
    "gene_name",
    "hpa_gene_name",
    "ensembl_id",
    "rna_tissue_specificity",
    "rna_tissue_distribution",
    "rna_tissue_specificity_score",
    "rna_tissue_specific_ntpm",
    "rna_single_cell_specificity",
    "rna_single_cell_distribution",
    "rna_single_cell_specificity_score",
    "rna_single_cell_specific_ntpm",
    "hpa_status",
]
thread_local = threading.local()


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Consulta a especificidade tecidual e de célula única "
            "dos genes no Human Protein Atlas usando gene_name, "
            "com consultas paralelas."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="CSV de entrada emitido pelo Nextflow.",
    )

    parser.add_argument(
        "--output",
        default="candidatos_HPA_tissue_specificity.csv",
        help="Nome do CSV de saída.",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Número de consultas simultâneas (padrão: 8).",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=20,
        help="Timeout de cada requisição em segundos.",
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=2,
        help="Número de tentativas adicionais em caso de falha.",
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=0,
        help="Espera opcional após cada consulta, em segundos.",
    )

    return parser.parse_args()


# HTTP SESSION
def get_session():
    """Retorna uma sessão HTTP reutilizável pela thread atual."""

    if not hasattr(thread_local, "session"):
        session = requests.Session()

        session.headers.update({
            "User-Agent": "tumor-antigen-pipeline/1.0",
        })

        thread_local.session = session

    return thread_local.session


# HPA QUERY
def query_hpa(gene_name, timeout=20, retries=2, delay=0):
    params = {
        "search": gene_name,
        "format": "json",
        "columns": HPA_COLUMNS,
        "compress": "no",
    }

    session = get_session()

    for attempt in range(retries + 1):
        try:
            response = session.get(
                HPA_API_URL,
                params=params,
                timeout=timeout,
            )

            if response.status_code == 429 or response.status_code >= 500:
                response.raise_for_status()

            response.raise_for_status()

            data = response.json()

            if isinstance(data, dict):
                data = [data]

            if not isinstance(data, list):
                return gene_name, None, "invalid_response"

            exact_matches = [
                item
                for item in data
                if isinstance(item, dict)
                and str(item.get("Gene", "")).strip().upper()
                == gene_name.strip().upper()
            ]

            if len(exact_matches) == 1:
                result = exact_matches[0]

                if delay > 0:
                    time.sleep(delay)

                return gene_name, result, "found"

            if len(exact_matches) > 1:
                return gene_name, None, "multiple_exact_matches"

            if data:
                return gene_name, None, "no_exact_match"

            return gene_name, None, "not_found"

        except (requests.RequestException, ValueError, json.JSONDecodeError) as exc:
            if attempt >= retries:
                return gene_name, None, f"error: {str(exc)[:200]}"
            time.sleep(min(2 ** attempt, 5))

    return gene_name, None, "request_failed"


# RESULT FORMATTING
def format_result(gene_name, hpa_data, status):
    result = {
        column: None
        for column in OUTPUT_COLUMNS
    }

    result["gene_name"] = gene_name
    result["hpa_status"] = status

    if hpa_data is None:
        return result

    result.update({
        "hpa_gene_name": hpa_data.get("Gene"),
        "ensembl_id": hpa_data.get("Ensembl"),

        "rna_tissue_specificity": hpa_data.get(
            "RNA tissue specificity"
        ),

        "rna_tissue_distribution": hpa_data.get(
            "RNA tissue distribution"
        ),

        "rna_tissue_specificity_score": hpa_data.get(
            "RNA tissue specificity score"
        ),

        "rna_tissue_specific_ntpm": hpa_data.get(
            "RNA tissue specific nTPM"
        ),

        "rna_single_cell_specificity": hpa_data.get(
            "RNA single cell type specificity"
        ),

        "rna_single_cell_distribution": hpa_data.get(
            "RNA single cell type distribution"
        ),

        "rna_single_cell_specificity_score": hpa_data.get(
            "RNA single cell type specificity score"
        ),

        "rna_single_cell_specific_ntpm": hpa_data.get(
            "RNA single cell type specific nCPM"
        ),
    })

    return result


# SAVE RESULTS
def save_results(results_by_gene, genes, output_path):
    ordered_results = [
        results_by_gene[gene]
        for gene in genes
        if gene in results_by_gene
    ]

    pd.DataFrame(
        ordered_results,
        columns=OUTPUT_COLUMNS,
    ).to_csv(output_path, index=False)


def main():
    args = parse_args()

    if args.workers < 1:
        raise ValueError("--workers precisa ser maior ou igual a 1.")

    if args.retries < 0:
        raise ValueError("--retries não pode ser negativo.")

    input_path = Path(args.input)
    output_path = Path(args.output)


    df = pd.read_csv(input_path)

    if "gene_name" not in df.columns:
        raise ValueError(
            "Fail: gene_name\n"
        )

    df["gene_name"] = df["gene_name"].astype("string").str.strip()

    genes = (
        df.loc[
            df["gene_name"].notna()
            & df["gene_name"].ne(""),
            "gene_name",
        ]
        .drop_duplicates()
        .tolist()
    )

    if not genes:
        raise ValueError(
            "No valid gene found in the column gene_name."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    results_by_gene = {}
    status_counter = Counter()
    completed = 0
    total = len(genes)

    with ThreadPoolExecutor(
        max_workers=args.workers
    ) as executor:

        futures = {
            executor.submit(
                query_hpa,
                gene_name,
                args.timeout,
                args.retries,
                args.delay,
            ): gene_name
            for gene_name in genes
        }

        try:
            for future in as_completed(futures):
                requested_gene = futures[future]

                try:
                    gene_name, hpa_data, status = future.result()

                except Exception as exc:
                    gene_name = requested_gene
                    hpa_data = None
                    status = f"error: {str(exc)[:200]}"

                result = format_result(
                    gene_name,
                    hpa_data,
                    status,
                )

                results_by_gene[gene_name] = result
                status_counter[status] += 1
                completed += 1

                save_results(
                    results_by_gene,
                    genes,
                    output_path,
                )

        except KeyboardInterrupt:
            for future in futures:
                future.cancel()

            save_results(
                results_by_gene,
                genes,
                output_path,
            )

            raise

    # SUMMARY
    results_df = pd.DataFrame(
        [results_by_gene[gene] for gene in genes],
        columns=OUTPUT_COLUMNS,
    )

    found_df = results_df[
        results_df["hpa_status"].eq("found")
    ]

    if not found_df.empty:
        specificity_counts = (
            found_df["rna_tissue_specificity"]
            .fillna("Not available")
            .value_counts()
        )


if __name__ == "__main__":
    main()
