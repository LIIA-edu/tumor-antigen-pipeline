#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze candidate gene heterogeneity in malignant cells."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--gene-map", required=True)
    parser.add_argument("--output_dir", default=".")
    parser.add_argument("--disease-column", default="disease")
    parser.add_argument("--cell-type-column", default="cell_type")
    parser.add_argument("--malignant-label", default="malignant cell")
    parser.add_argument("--gene-column", default="gene_name")
    parser.add_argument("--n-modules", type=int, default=10)

    return parser.parse_args()


def to_dense(matrix):
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    candidates_df = pd.read_csv(args.candidates)

    with open(args.gene_map, "r", encoding="utf-8") as handle:
        gene_map = json.load(handle)

    if args.gene_column not in candidates_df.columns:
        raise ValueError(
            f"Candidate column '{args.gene_column}' not found. "
            f"Available columns: {list(candidates_df.columns)}"
        )

    for column in (args.disease_column, args.cell_type_column):
        if column not in adata.obs.columns:
            raise ValueError(
                f"Column '{column}' not found in AnnData.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    # Keep malignant cells
    cell_types = (
        adata.obs[args.cell_type_column]
        .astype("string")
        .str.strip()
        .str.lower()
    )

    malignant_mask = (
        cell_types == args.malignant_label.strip().lower()
    ).to_numpy()

    if not malignant_mask.any():
        raise ValueError(
            f"No cells found with {args.cell_type_column}="
            f"'{args.malignant_label}'."
        )

    malignant = adata[malignant_mask].copy()

    # Map feature IDs to gene symbols
    mapped_genes = []

    for feature_id in malignant.var_names.astype(str):
        gene = gene_map.get(feature_id, feature_id)
        if gene.startswith("ENSG"):
            gene = gene.split(".")[0]
        mapped_genes.append(gene)

    # Match candidate genes
    candidate_genes = (
        candidates_df[args.gene_column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    candidate_genes = {
        gene for gene in candidate_genes if gene
    }

    gene_indices = {}

    for index, gene in enumerate(mapped_genes):
        if gene in candidate_genes:
            gene_indices.setdefault(gene, []).append(index)

    matched_genes = sorted(gene_indices)

    if not matched_genes:
        raise ValueError(
            "No candidate genes matched the AnnData feature names. "
            "Check the gene mapping and candidate gene column."
        )

    # Calculate expression for each candidate gene
    X = malignant.X

    diseases = (
        malignant.obs[args.disease_column]
        .astype("string")
        .fillna("Unknown")
        .astype(str)
        .to_numpy()
    )

    gene_expression = {}

    for gene in matched_genes:
        indices = gene_indices[gene]
        values = to_dense(X[:, indices])

        if values.ndim == 1:
            values = values.reshape(-1, 1)

        # Combine features that map to the same gene.
        if values.shape[1] > 1:
            values = values.sum(axis=1)
        else:
            values = values[:, 0]

        gene_expression[gene] = np.asarray(
            values, dtype=float
        ).reshape(-1)

    
    # Calculate expression by disease
    expression_records = []

    for disease in sorted(set(diseases)):
        disease_mask = diseases == disease
        n_cells = int(disease_mask.sum())

        for gene in matched_genes:
            values = gene_expression[gene][disease_mask]

            expression_records.append({
                "gene": gene,
                "disease": disease,
                "mean_expression": float(np.mean(values)),
                "pct_expressing": float(
                    np.mean(values > 0) * 100
                ),
                "n_malignant_cells": n_cells,
            })

    expression_df = pd.DataFrame(expression_records)

    # Generate mean expression matrix
    # Rows = genes; columns = diseases
    expression_matrix = expression_df.pivot(
        index="gene",
        columns="disease",
        values="mean_expression",
    ).fillna(0.0)

    expression_matrix = expression_matrix.reindex(
        matched_genes,
        fill_value=0.0,
    )

    mean_expression_file = (
        output_dir /
        "candidate_genes_mean_expression_by_disease.csv"
    )

    expression_matrix.to_csv(
        mean_expression_file,
        index=True,
        index_label="gene",
    )
    
    # Generate percentage expression matrix
    # Rows = genes; columns = diseases
    # Values = percentage of malignant cells expressing each gene
    pct_expression_matrix = expression_df.pivot(
        index="gene",
        columns="disease",
        values="pct_expressing",
    ).fillna(0.0)

    pct_expression_matrix = pct_expression_matrix.reindex(
        matched_genes,
        fill_value=0.0,
    )

    pct_expression_file = (
        output_dir /
        "candidate_genes_pct_expression_by_disease.csv"
    )

    pct_expression_matrix.to_csv(
        pct_expression_file,
        index=True,
        index_label="gene",
    )
    
    # Assign genes to expression modules
    n_genes = len(expression_matrix)

    if n_genes == 1:
        assignments = pd.DataFrame({
            "gene": expression_matrix.index,
            "module": ["Module_1"],
        })
    else:
        n_modules = max(1, min(args.n_modules, n_genes))

        scaled_expression = StandardScaler().fit_transform(
            expression_matrix.to_numpy()
        )

        model = KMeans(
            n_clusters=n_modules,
            random_state=42,
            n_init=10,
        )

        labels = model.fit_predict(scaled_expression)

        assignments = pd.DataFrame({
            "gene": expression_matrix.index,
            "module": [
                f"Module_{label + 1}" for label in labels
            ],
        })

    # Generate candidate_gene_modules.csv
    candidate_gene_modules = expression_df.merge(
        assignments,
        on="gene",
        how="inner",
    )

    candidate_gene_modules = candidate_gene_modules[
        [
            "module",
            "gene",
            "disease",
            "mean_expression",
            "pct_expressing",
            "n_malignant_cells",
        ]
    ]

    modules_file = output_dir / "candidate_gene_modules.csv"
    candidate_gene_modules.to_csv(
        modules_file,
        index=False,
    )

    # Report
    for file_path in (
        modules_file,
        mean_expression_file,
        pct_expression_file,
    ):
        print(f"  - {file_path}", flush=True)


if __name__ == "__main__":
    main()
