#!/usr/bin/env python3

import argparse
import pandas as pd
from antigens_type import ANTIGENS_TYPE

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
args = parser.parse_args()
df = pd.read_csv(args.input)


# CLASSIFY ANTIGEN TYPES
for type_antigen, list_antigens in ANTIGENS_TYPE.items():

    mask = df["gene_name"].isin(list_antigens)
    df.loc[mask, "type"] = type_antigen


# SAVE CLASSIFIED SPECIFICITY
df.to_csv(
    "classified_specificity_labeled.csv",
    index=False
)


# SAVE MERGED SUMMED GENES
merged_summed_genes = df[
    df["type"].notna()
].copy()

merged_summed_genes = merged_summed_genes[
    ["gene_name", "type"]
].drop_duplicates()

merged_summed_genes.to_csv(
    "merged_summed_genes.csv",
    index=False
)
