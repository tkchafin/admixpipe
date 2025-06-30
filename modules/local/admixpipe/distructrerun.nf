process DISTRUCT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/mussmann/admixpipe:3.2'

    input:
    tuple val(meta), path(pfiles)
    tuple val(meta2), path(qfiles)
    tuple val(meta3), path(pops)
    tuple val(meta4), path(inds)
    tuple val(meta5), path(logs)
    tuple val(meta6), path(clumpak)

    output:
    tuple val(meta), path("MajorClusterRuns.txt"), emit: major_clusters
    tuple val(meta), path("*/best_results"), emit: best_results
    tuple val(meta), path("cv_file.MajClust.txt"), emit: cv
    tuple val(meta), path("loglikelihood_file.MajClust.txt"), emit: loglik
    tuple val(meta), path("cvRuns.json"), emit: cvruns_json
    tuple val(meta), path("qfilePaths.json"), emit: qfilepaths_json
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml",  emit: versions

    script:
    def args   = task.ext.args ?: ''
    def maxk   = params.maxk ?: 3

    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/clumpak:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    distructRerun.py \\
    -a ./ \\
    -d ${clumpak} \\
    -k 1 \\
    -K ${maxk} \\
    -r \\
    ${args}

    for f in */best_results/*.ps; do
        ps2pdf \$f;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2
        distruct: 1.1
        ghostscript: 9.50
    END_VERSIONS
    """
}
