process PLOT_EVALADMIX {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(qfilepaths)
        tuple val(meta2), path(fam)
        tuple val(meta3), path(corres)
        tuple val(meta4), path(bestk_file)

    output:
        path("evaladmix_allk_mqc.html"), emit: allk_html
        path("evaladmix_bestk_mqc.html"), emit: bestk_html
        path("versions.yml")   , emit: versions

    script:
    """
    bestk=\$(cat ${bestk_file})

    plot_evaladmix.py \\
            --prefix "${meta.id}" \\
            --workdir . \\
            --qfilePaths ${qfilepaths} \\
            --allk_html evaladmix_allk_mqc.html \\
            --template_allk ${baseDir}/assets/multiqc_evaladmix_allk.html \\
            --best_k \$bestk \\
            --bestk_html evaladmix_bestk_mqc.html \\
            --template_bestk ${baseDir}/assets/multiqc_evaladmix_bestk.html

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
