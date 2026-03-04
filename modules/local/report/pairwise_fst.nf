process PLOT_PAIRWISE_FST {
    label 'process_single'
    tag "$meta.id"

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(snpio_report_data)

    output:
        path("pairwise_fst_mqc.html"), emit: plot_html
        path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    echo "🔍 Finding input files..."
    fst=\$(find -L ${snpio_report_data} -type f -name 'wc_fst_permutation_observed.txt' | head -n1)
    pvals=\$(find -L ${snpio_report_data} -type f -name 'wc_fst_permutation_pvalues.txt' | head -n1)

    echo "📊 Plotting..."
    plot_pairwise_fst.py \\
        --fst \$fst \\
        --pvals \$pvals \\
        --header ${baseDir}/assets/multiqc_pairwise_fst.html \\
        --output pairwise_fst_mqc.html

    pandas_version=\$(python3 -c 'import pandas; print(pandas.__version__)')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \${pandas_version}
    END_VERSIONS
    """
}
