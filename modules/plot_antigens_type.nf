
process PLOT_ANTIGENS_TYPE {

    tag "${specificity_csv.simpleName}"

    publishDir "${params.outdir}/antigens_type",
        mode: 'copy',
        overwrite: true

    input:
    path specificity_csv

    output:
    path "specificity_groups_by_type_heatmap.svg", emit: heatmap
    path "merged_summed_genes.csv", emit: merged
    path "cta_malignant_vs_others_venn.svg", emit: venn
    path "cta_enriched_categories_barplot.svg", emit: barplot

    script:
    """
    Rscript "${projectDir}/scripts/plot_antigens_type.R" \
        --input "${specificity_csv}" \
        --output specificity_groups_by_type_heatmap.svg \
        --merged_output merged_summed_genes.csv \
        --venn_output cta_malignant_vs_others_venn.svg \
        --barplot_output cta_enriched_categories_barplot.svg
    """
}
