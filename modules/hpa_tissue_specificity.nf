process HPA_TISSUE_SPECIFICITY {

    tag "HPA tissue specificity"

    publishDir "${params.outdir ?: 'results'}/hpa_tissue_specificity",
        mode: 'copy',
        overwrite: true

    input:
    path candidates_csv

    output:
    path 'candidatos_HPA_tissue_specificity.csv',
        emit: tissue_specificity

    script:
    """
    python3 "${projectDir}/scripts/hpa_tissue_specificity.py" \
        --input "${candidates_csv}" \
        --output candidatos_HPA_tissue_specificity.csv
    """
}
