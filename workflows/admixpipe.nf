/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_admixpipe_pipeline'

include { SNPIO_FILTER } from '../modules/local/snpio/filter.nf'
include { RUN_ADMIXPIPE } from '../subworkflows/local/run_admixpipe.nf'
include { GENERATE_REPORT } from '../subworkflows/local/generate_report.nf'
include { CUSTOMIZE_REPORT } from '../modules/local/report/customize_report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ADMIXPIPE {

    take:
    ch_vcf     // [meta, vcf]
    ch_tbi     // [meta, tbi]
    ch_popmap  // [meta, popmap]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // VCF pre-processing
    //
    // This step removes individuals with a large amount of missing data,
    // flanking variants, low variation, etc...
    // and also generates SNPio missingness reports
    SNPIO_FILTER(
        ch_vcf,
        ch_tbi,
        ch_popmap
    )
    ch_versions = ch_versions.mix(SNPIO_FILTER.out.versions)
    ch_filtered_vcf = SNPIO_FILTER.out.filtered_vcf.map { meta, file -> tuple(meta + [id: "${meta.id}_filtered"], file) }
    ch_filtered_tbi = SNPIO_FILTER.out.filtered_tbi.map { meta, file -> tuple(meta + [id: "${meta.id}_filtered"], file) }
    ch_snpio_output = SNPIO_FILTER.out.snpio_output.map { meta, dir -> tuple(meta + [id: "${meta.id}_filtered"], dir) }

    //
    // Run admixture pipeline on filtered dataset
    //
    RUN_ADMIXPIPE(
        ch_filtered_vcf,
        ch_popmap
    )
    ch_versions = ch_versions.mix(RUN_ADMIXPIPE.out.versions)


    //
    // Generate figures for the report
    //
    GENERATE_REPORT(
        ch_vcf,
        ch_tbi,
        ch_filtered_vcf,
        ch_filtered_tbi,
        RUN_ADMIXPIPE.out.cv_file,
        RUN_ADMIXPIPE.out.evanno,
        RUN_ADMIXPIPE.out.bestK_file,
        RUN_ADMIXPIPE.out.best_results,
        ch_snpio_output,
        RUN_ADMIXPIPE.out.bestK_clumpp,
        RUN_ADMIXPIPE.out.inds,
        RUN_ADMIXPIPE.out.pops,
    )
    ch_versions = ch_versions.mix( GENERATE_REPORT.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( GENERATE_REPORT.out.mqc_files )


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // Customize report
    CUSTOMIZE_REPORT( MULTIQC.out.report )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
