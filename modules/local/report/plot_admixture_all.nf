process PLOT_ADMIXTURE_ALL {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(best_results)
        tuple val(meta2), path(inds)
        tuple val(meta3), path(pops)
        tuple val(meta4), path(bestk_file)
    output:
        path("admixture_allk_mqc.html"), emit: admixture_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    bestk=\$(cat ${bestk_file})

    plot_admixture_all.py \\
        --indir ${best_results} \\
        --inds ${inds} \\
        --pops ${pops} \\
        --template ${baseDir}/assets/multiqc_admixture_allk.html \\
        --out "admixture_allk_mqc.html" \\
        --bestk \$bestk \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
