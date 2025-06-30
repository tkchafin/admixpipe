process PLOT_CV {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(cv_file)

    output:
        path("cvplot_mqc.html"), emit: cv_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_cv.py \\
        ${cv_file} \\
        --template ${baseDir}/assets/multiqc_cv.html \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
