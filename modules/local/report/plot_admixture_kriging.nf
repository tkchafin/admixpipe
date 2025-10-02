process PLOT_ADMIXTURE_KRIGING {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/pykrige:1.0"

    input:
    tuple val(meta),  path(clumppfile)    // .Q file
    tuple val(meta2), path(inds)          // one sample ID per line
    tuple val(meta3), path(pops)          // one pop/site ID per line
    tuple val(meta4), path(site_coords)   // TSV: ID,Lat,Lon
    tuple val(meta5), path(geo_data)      // optional dir with config.json

    output:
    path("admixture_kriging_continuous.html"),         emit: html_continuous
    path("admixture_kriging_continuous_merged.html"),  emit: html_continuous_merged
    path("admixture_kriging_discrete.html"),           emit: html_discrete
    path("admixture_kriging_simpson.html"),            emit: html_simpson
    path("admixture_kriging_continuous.tif"),          emit: tif_continuous
    path("admixture_kriging_discrete.tif"),            emit: tif_discrete
    path("admixture_kriging_simpson.tif"),             emit: tif_simpson
    path("admixture_kriging_continuous_merged.tif"),  emit: tif_merged
    path("versions.yml"),                               emit: versions

    script:
    def args         = task.ext.args ?: ''
    def geo_data_arg = geo_data ? "--geo_data_json ${geo_data}/config.json" : ''
    """
    kriging_singlek.py \\
        --qmat         ${clumppfile} \\
        --inds         ${inds} \\
        --pops         ${pops} \\
        --site_coords  ${site_coords} \\
        --out_prefix   admixture \\
        --template_cont   ${baseDir}/assets/multiqc_kriging_continuous.html \\
        --template_disc   ${baseDir}/assets/multiqc_kriging_discrete.html \\
        --template_div    ${baseDir}/assets/multiqc_kriging_simpson.html \\
        --template_merged ${baseDir}/assets/multiqc_kriging_continuous_merged.html \\
        --grid_nx 500 --grid_ny 500 \\
        ${geo_data_arg} \\
        ${args}

    python_version=\$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
    numpy_version=\$(python3 -c 'import numpy as np; print(np.__version__)')
    pykrige_version=\$(python3 -c 'import pykrige; print(pykrige.__version__)')
    geopandas_version=\$(python3 -c 'import geopandas; print(geopandas.__version__)')
    rasterio_version=\$(python3 -c 'import rasterio; print(rasterio.__version__)')
    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:    \${python_version}
        numpy:     \${numpy_version}
        pykrige:   \${pykrige_version}
        geopandas: \${geopandas_version}
        rasterio:  \${rasterio_version}
        plotly:    \${plotly_version}
    END_VERSIONS
    """
}
