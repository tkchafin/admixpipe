process EVANNO {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
        tuple val(meta), path(ll_file)

    output:
        tuple val(meta), path('evanno_metrics.tsv')  , emit: metrics

    script:
    """
    evanno.py \\
        ${ll_file} \\
        -o "evanno_metrics.tsv"
    """
}
