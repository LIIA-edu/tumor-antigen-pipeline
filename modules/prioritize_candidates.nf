process PRIORITIZE_CANDIDATES {

    tag "Candidate prioritization"

    publishDir "${params.outdir ?: 'results'}/prioritized_candidates",
        mode: 'copy',
        overwrite: true

    cpus 1
    memory '2 GB'
    time '1h'

    input:
    path candidates_csv
    path hpa_csv
    path epitopes_csv
    path binary_matrix
    path merged_summed_genes

    output:
    path 'prioritized_candidates.csv',
        emit: prioritized_candidates

    script:
    """
    python3 "${projectDir}/scripts/prioritize_candidates.py" \
        --candidates "${candidates_csv}" \
        --hpa "${hpa_csv}" \
        --epitopes "${epitopes_csv}" \
        --binary-matrix "${binary_matrix}" \
        --merged-summed-genes "${merged_summed_genes}" \
        --output prioritized_candidates.csv
    """
}