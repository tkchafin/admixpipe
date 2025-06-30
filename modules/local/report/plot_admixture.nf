process PLOT_ADMIXTURE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(clumppfile), path(inds), path(pops)
    output:
        path("admixture_${meta.id}_mqc.html"), emit: admixture_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_admixture.py \\
        --clumpp ${clumppfile} \\
        --inds ${inds} \\
        --pops ${pops} \\
        --template ${baseDir}/assets/multiqc_admixture_${meta.id}.html \\
        --out "admixture_${meta.id}_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
