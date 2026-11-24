process PLOT_EPITOPES {

    tag "epitope_plots"

    publishDir "${params.outdir}/epitope_analysis",
        mode: 'copy',
        overwrite: true

    input:
    path epitopes
    path candidates

    output:
    path "plot_A_genes_with_known_epitope_associations.svg"
    path "plot_B_cancer_associations.svg"
    path "genes_by_epitope_count.csv"

    script:
    """
    Rscript "${projectDir}/scripts/plot_epitopes.R" \
        "${epitopes}" \
        "${candidates}" \
        .
    """
}