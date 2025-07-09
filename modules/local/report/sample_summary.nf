process SAMPLE_SUMMARY {
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(snpio_report_data)
        tuple val(meta2), path(inds_pre)
        tuple val(meta3), path(inds_post)

    output:
        path("sample_summary_mqc.json"), emit: summary_txt
        path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    echo "🔍 Finding input files..."
    miss1=\$(find -L ${snpio_report_data} -type f -name 'multiqc_individual_missingness_table.txt' | head -n1)
    miss2=\$(find -L ${snpio_report_data} -type f -name 'multiqc_individual_missingness-2_table.txt' | head -n1)

    echo "📊 Generating sample summary..."
    sample_summary.py \\
        --inds-pre ${inds_pre} \\
        --inds-post ${inds_post} \\
        --miss-pre \$miss1 \\
        --miss-post \$miss2 \\
        --header ${baseDir}/assets/multiqc_sample_stats.html \\
        --output sample_summary_mqc.json

    pandas_version=\$(python3 -c 'import pandas; print(pandas.__version__)')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \${pandas_version}
    END_VERSIONS
    """
}
