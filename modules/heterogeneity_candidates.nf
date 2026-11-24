process HETEROGENEITY_CANDIDATES {

    tag "heterogeneity_candidates"

    publishDir "${params.outdir}/heterogeneity_candidates",
        mode: 'copy',
        overwrite: true

    input:
    path h5ad
    path candidates
    path gene_map

    output:
    path "candidate_gene_modules.csv", emit: modules_by_disease
    path "candidate_genes_mean_expression_by_disease.csv", emit: mean_expression_by_disease
    path "candidate_genes_pct_expression_by_disease.csv", emit: pct_expression_by_disease

    script:
    """
    python "${projectDir}/scripts/heterogeneity_candidates.py" \
        --input "${h5ad}" \
        --candidates "${candidates}" \
        --gene-map "${gene_map}" \
        --output_dir "." \
        --disease-column disease \
        --cell-type-column cell_type \
        --malignant-label "malignant cell" \
        --gene-column gene_name
    """
}