process EVALADMIX {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/mussmann/admixpipe:3.2'

    input:
    tuple val(meta), path(ped)
    tuple val(meta2), path(map)
    tuple val(meta3), path(pfiles)
    tuple val(meta4), path(qfiles)
    tuple val(meta5), path(qfiles_json)
    tuple val(meta6), path(popmap)
    tuple val(meta7), path(clumpak_output)
    tuple val(meta8), path(major_clusters)
    tuple val(meta9), path(cvruns_json)
    tuple val(meta10), path(qfilepaths_json)

    output:
    tuple val(meta), path("*corres"), emit: corres
    tuple val(meta), path("[1-9]*.png"), emit: majorclust_png
    tuple val(meta), path("*MinClust*.png"), optional:true, emit: minorclust_png
    tuple val(meta), path("${meta.id}*.png"), emit: reps_png
    tuple val(meta), path("${meta.id}*.fam"), emit: fam
    path "versions.yml",  emit: versions

    script:
    def args   = task.ext.args ?: ''

    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/clumpak:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    runEvalAdmix.py \\
    -p ${meta.id} \\
    -k 1 \\
    -K ${params.maxk} \\
    -m ${popmap} \\
    -n ${task.cpus} \\
    ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2
    END_VERSIONS
    """
}
