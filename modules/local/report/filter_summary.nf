process FILTER_SUMMARY {
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
        tuple val(meta), path(snpio_pre)

    output:
        path("sankey_mqc.html"), emit: sankey_html

    script:
    def args = task.ext.args ?: ''

    """
    html=\$(find -L ${snpio_pre} -type f -name 'filtering_results_sankey*.html' | head -n1)

    cat ${baseDir}/assets/multiqc_sankey.html > sankey_mqc.html
    cat \$html >> sankey_mqc.html
    """
}
