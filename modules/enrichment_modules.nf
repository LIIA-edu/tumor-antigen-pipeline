
process ENRICHMENT_MODULES {

    tag "Functional enrichment"

    publishDir "${params.outdir ?: 'results'}/enrichment",
        mode: 'copy',
        overwrite: true

    cpus 1
    memory '2 GB'
    time '2h'

    input:
    path input_csv

    output:
    path 'enrichment_results',
        emit: enrichment_results

    script:
    """
    python3 "${projectDir}/scripts/enrichment_modules.py" \
        --input "${input_csv}" \
        --outdir enrichment_results
    """
}
