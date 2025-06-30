process SAMPLE_SUMMARY {
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(inds_pre)
        tuple val(meta2), path(inds_post)
        tuple val(meta3), path(snpio_pre)
        tuple val(meta4), path(snpio_post)

    output:
        path("sample_summary_mqc.json"), emit: summary_txt
        path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    echo "🔍 Finding input files..."
    miss1=\$(find -L ${snpio_pre} -type f -name 'individual_missingness.csv' | head -n1)
    miss2=\$(find -L ${snpio_post} -type f -name 'individual_missingness.csv' | head -n1)
    het1=\$(find -L ${snpio_pre} -type f -name 'pop_individ_locus_missingness.csv' | head -n1)
    het2=\$(find -L ${snpio_post} -type f -name 'pop_individ_locus_missingness.csv' | head -n1)

    echo "📊 Generating sample summary..."
    sample_summary.py \\
        --inds-pre ${inds_pre} \\
        --inds-post ${inds_post} \\
        --miss-pre \$miss1 \\
        --miss-post \$miss2 \\
        --het-pre \$het1 \\
        --het-post \$het2 \\
        --header ${baseDir}/assets/multiqc_sample_stats.html \\
        --output sample_summary_mqc.json

    pandas_version=\$(python3 -c 'import pandas; print(pandas.__version__)')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \${pandas_version}
    END_VERSIONS
    """
}
