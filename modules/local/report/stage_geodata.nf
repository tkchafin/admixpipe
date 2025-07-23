process STAGE_GEODATA_LAYERS {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(json)
        tuple val(meta), path(dir)
    output:
        tuple val(meta), path("geo_data_files"), emit: geo_data_dir

    script:
    def args   = task.ext.args ?: ''
    """
    stage_geodata.py \\
        --json   ${json} \\
        --outdir geo_data_files \\
        --out    geo_data_files/config.json \\
        --base_path . \\
        ${args}
    """
}
