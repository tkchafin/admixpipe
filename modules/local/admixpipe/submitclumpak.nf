process CLUMPAK {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/mussmann/admixpipe:3.2'

    input:
    tuple val(meta), path(results)
    tuple val(meta2), path(inds)
    tuple val(meta3), path(pops)

    output:
    tuple val(meta), path("clumpakOutput"),   emit: output
    path "versions.yml",  emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/clumpak:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    submitClumpak.py \\
        -r ${results} \\
        -p ${meta.id} \\
        -M \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2
        CLUMPAK: 1.1
    END_VERSIONS
    """
}
