process PLOT_PCA {
    label 'process_single'
    tag "$meta.id"

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(snpio_report_data)

    output:
        path("snpio_pca_mqc.html"), emit: plot_html
        path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    echo "🔍 Finding input files..."
    pca=\$(find -L ${snpio_report_data} -type f -name 'table*.txt' | head -n1)

    echo "📊 Plotting..."
    plot_pca.py \\
        --input \$pca \\
        --header ${baseDir}/assets/multiqc_pca.html \\
        --output snpio_pca_mqc.html

    pandas_version=\$(python3 -c 'import pandas; print(pandas.__version__)')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \${pandas_version}
    END_VERSIONS
    """
}
