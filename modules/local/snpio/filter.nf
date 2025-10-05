process SNPIO_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/btmartin721/snpio:1.6.10'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(popmap)

    output:
    tuple val(meta), path("${meta.id}.filter.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.filter.vcf.gz.tbi"), emit: filtered_tbi
    tuple val(meta), path("*_output"), emit: snpio_output
    tuple val(meta), path("multiqc/multiqc_report.html"), emit: multiqc_report
    tuple val(meta), path("multiqc/multiqc_report_data"), emit: multiqc_report_data
    path "versions.yml",     emit: versions

    script:
    def args   = task.ext.args ?: ''

    """
    snpio_filter.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --ind_cov ${params.ind_cov} \\
        --min_maf ${params.min_maf} \\
        --snp_cov ${params.snp_cov} \\
        --ind_cov ${params.ind_cov} \\
        --pop_cov ${params.pop_cov} \\
        --flank_dist ${params.thin_dist} \\
        --jobs ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SNPio: 1.6.10
    END_VERSIONS
    """
}
