process PLOT_HETEROGENEITY_CANDIDATES {

    tag "plot_heterogeneity_candidates"

    publishDir "${params.outdir ?: 'results'}/heterogeneity_candidates",
        mode: 'copy',
        overwrite: true

    input:
    path candidate_gene_modules
    path mean_expression_by_disease
    path pct_expression_by_disease
    path epitopes_file

    output:
    path "*.png", optional: true
    path "*.pdf", optional: true
    path "*.svg", optional: true
    path "*.csv", optional: true

    script:
    """
    Rscript "${projectDir}/scripts/plot_heterogeneity_candidates.R" \
        --modules "${candidate_gene_modules}" \
        --mean-expression "${mean_expression_by_disease}" \
        --pct-expression "${pct_expression_by_disease}" \
        --epitopes "${epitopes_file}"
    """
}