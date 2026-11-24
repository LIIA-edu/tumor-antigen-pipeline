
process DATA_CATEGORIZATION {

    tag "${h5ad.simpleName}"

    publishDir "${params.outdir}/data_categorization",
        mode: 'copy',
        overwrite: true

    input:
    path h5ad

    output:
    path "Figure1C_CellsPerTumorType.svg"
    path "Figure1C_CellsPerTumorType.png"
    path "Figure1D_CellTypeComposition.svg"
    path "Figure1D_CellTypeComposition.png"
    path "Figure1E_MalignantVsNonMalignant_UMAP.svg"
    path "Figure1E_MalignantVsNonMalignant_UMAP.png"
    path "Figure1F_Disease_UMAP.png"

    script:
    """
    Rscript ${projectDir}/scripts/plot_data_categorization.R \
        --input-h5ad ${h5ad} \
        --output-dir .
    """
}
