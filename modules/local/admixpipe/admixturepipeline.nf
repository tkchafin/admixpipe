process ADMIXTUREPIPELINE {
    tag "$meta.id"
    label 'process_large'

    container 'docker.io/mussmann/admixpipe:3.2'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)

    output:
    tuple val(meta), path("results.zip"),      emit: results
    tuple val(meta), path("${meta.id}*.stdout"),         emit: logs
    tuple val(meta), path("${meta.id}*.Q"),              emit: qfiles
    tuple val(meta), path("${meta.id}*.P"),              emit: pfiles
    tuple val(meta), path("${meta.id}_pops.txt"),        emit: pops
    tuple val(meta), path("${meta.id}_inds.txt"),        emit: inds
    tuple val(meta), path("${meta.id}.map"),             emit: map
    tuple val(meta), path("${meta.id}.ped"),             emit: ped
    tuple val(meta), path("${meta.id}.qfiles.json"),     emit: qfiles_json
    path "versions.yml",     emit: versions

    script:
    def args   = task.ext.args ?: ''
    def maxk   = params.maxk ?: 3

    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    admixturePipeline.py \\
        -m ${popmap} \\
        -v ${vcf} \\
        -n ${task.cpus} \\
        -k 1 \\
        -K ${maxk} \\
        -C 1.0 \\
        -S 0.0 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2
        VCFtools: 0.1.16
        PLINK: 20220402
        Admixture: 1.30
    END_VERSIONS
    """
}
