process PLOT_EVANNO {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(evanno_file)
        tuple val(meta2), path(bestk_file)

    output:
        path("evanno_mqc.html"), emit: evanno_html
        path("versions.yml"),     emit: versions

    script:
    """
    bestk=\$(cat ${bestk_file})

    plot_evanno.py \\
        ${evanno_file} \\
        --bestk \$bestk \\
        --template ${baseDir}/assets/multiqc_evanno.html \\
        -o evanno_mqc.html

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
