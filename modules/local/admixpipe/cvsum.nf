process CVSUM {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/mussmann/admixpipe:3.2'

    input:
    tuple val(meta), path(cv)
    tuple val(meta2), path(loglik)

    output:
    tuple val(meta), path("cv_file.MajClust.png"), emit: cv_plot
    tuple val(meta), path("loglikelihood_file.MajClust.png"), emit: loglik_plot
    tuple val(meta), path("cv_output.txt"), emit: cv_output
    tuple val(meta), path("ll_output.txt"), emit: ll_output
    path "versions.yml", emit: versions

    script:
    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/clumpak:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    cvSum.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2
    END_VERSIONS
    """
}
