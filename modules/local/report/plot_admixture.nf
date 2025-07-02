process PLOT_ADMIXTURE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(clumppfile)
        tuple val(meta2), path(inds)
        tuple val(meta3), path(pops)
    output:
        path("admixture_bestk_mqc.html"), emit: admixture_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_admixture.py \\
        --clumpp ${clumppfile} \\
        --inds ${inds} \\
        --pops ${pops} \\
        --template ${baseDir}/assets/multiqc_admixture_bestk.html \\
        --out "admixture_bestk_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
