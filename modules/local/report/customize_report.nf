process CUSTOMIZE_REPORT {
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
        path(in_report), stageAs: "in_report.html"

    output:
        path("multiqc_report.html"), emit: report

    script:
    """
    encoded_logo=\$(cat ${baseDir}/docs/images/logo.b64)

    echo "<h1 id=\\"page_title\\">
    <a href=\\"https://github.com/UARK-aCaMEL\\" target=\\"_blank\\">
        <img src=\\"data:image/png;base64,\${encoded_logo}\\" title=\\"Arkansas Conservation and Molecular Ecology Lab\\" class=\\"multiqc_logo\\">
    </a>
</h1>" > new_header.html

    awk '
        BEGIN {in_block=0}
        /<h1 id="page_title">/ {print "_HEADER_REPLACEMENT_"; in_block=1; next}
        /<\\/h1>/ && in_block {in_block=0; next}
        !in_block {print}
    ' in_report.html | sed -e "/_HEADER_REPLACEMENT_/ {
        r new_header.html
        d
    }" > multiqc_report.html

    rm new_header.html
    """
}
