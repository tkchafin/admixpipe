include { PLOT_CV } from '../../modules/local/report/plot_cv.nf'
include { PLOT_EVANNO } from '../../modules/local/report/plot_evanno.nf'
include { SAMPLE_SUMMARY } from '../../modules/local/report/sample_summary.nf'
include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'
include { PLOT_ADMIXTURE_ALL } from '../../modules/local/report/plot_admixture_all.nf'
include { FILTER_SUMMARY } from '../../modules/local/report/filter_summary.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_PRE } from '../../modules/local/bcftools_query.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_POST } from '../../modules/local/bcftools_query.nf'


workflow GENERATE_REPORT {
    take:
    vcf_pre
    tbi_pre
    vcf_post
    tbi_post
    cv_file
    evanno
    bestk_file
    best_results
    snpio_pre
    clumpp
    inds
    pops

    main:
    ch_versions = Channel.empty()
    ch_mqc_files = Channel.empty()

    //CV plot
    PLOT_CV( cv_file, bestk_file )
    ch_versions = ch_versions.mix( PLOT_CV.out.versions )
    ch_mqc_files = ch_mqc_files.mix( PLOT_CV.out.cv_html )

    //Evanno plot
    PLOT_EVANNO( evanno, bestk_file )
    ch_versions = ch_versions.mix( PLOT_EVANNO.out.versions )
    ch_mqc_files = ch_mqc_files.mix( PLOT_EVANNO.out.evanno_html )

    //Get individual lists from vcfs
    BCFTOOLS_QUERY_PRE( vcf_pre, tbi_pre )
    BCFTOOLS_QUERY_POST( vcf_post, tbi_post )
    ch_versions = ch_versions.mix( BCFTOOLS_QUERY_PRE.out.versions )

    // //SNPio summary
    // SAMPLE_SUMMARY(
    //     BCFTOOLS_QUERY_PRE.out.samples,
    //     BCFTOOLS_QUERY_POST.out.samples,
    //     snpio_pre
    // )
    // ch_mqc_files = ch_mqc_files.mix( SAMPLE_SUMMARY.out.summary_txt )
    // ch_versions = ch_versions.mix( SAMPLE_SUMMARY.out.versions )


    //Admixture barplots
    PLOT_ADMIXTURE(
        clumpp,
        inds,
        pops
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE.out.versions )

    //Admixture barplots -- all
    PLOT_ADMIXTURE_ALL(
        best_results,
        inds,
        pops,
        bestk_file
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_ALL.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE_ALL.out.versions )


    //SNPio plots
    FILTER_SUMMARY( snpio_pre )
    ch_mqc_files = ch_mqc_files.mix( FILTER_SUMMARY.out.sankey_html )

    emit:
    mqc_files    = ch_mqc_files
    versions     = ch_versions
}
