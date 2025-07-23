process PLOT_ADMIXTURE_SPATIAL {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/geopandas:1.0"

    input:
    tuple val(meta),   path(clumppfile)   // .Q file
    tuple val(meta2),  path(inds)         // one sample ID per line
    tuple val(meta3),  path(pops)         // one pop/site ID per line
    tuple val(meta4),  path(site_coords)  // TSV: ID,Lat,Lon
    tuple val(meta5), path(geo_data)

    output:
    path("admixture_spatial_mqc.html"), emit: plot_html
    path("admixture_spatial.tsv"),      emit: admixture_summary
    path("versions.yml"),               emit: versions

    script:
    def args   = task.ext.args ?: ''
    def geo_data_arg = geo_data ? "--geo_data_json ${geo_data}/config.json" : ''
    """
    plot_admixture_spatial.py \\
        --qmat       ${clumppfile} \\
        --inds       ${inds} \\
        --pops       ${pops} \\
        --site_coords ${site_coords} \\
        --template   ${baseDir}/assets/multiqc_admixture_spatial.html \\
        --out        admixture_spatial_mqc.html \\
        --table_out  admixture_spatial.tsv \\
        ${geo_data_arg} \\
        ${args}

    folium_version=\$(python3 -c 'import folium; print(folium.__version__)')
    geopandas_version=\$(python3 -c 'import geopandas; print(geopandas.__version__)')
    matplotlib_version=\$(python3 -c 'import matplotlib; print(matplotlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        folium:     \${folium_version}
        geopandas:  \${geopandas_version}
        matplotlib: \${matplotlib_version}
    END_VERSIONS
    """
}
