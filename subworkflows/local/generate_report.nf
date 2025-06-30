//
// Rank loci using Rosenberg et al. 2003 importance indices
//
include { PLOT_CV } from '../../modules/local/report/plot_cv.nf'
include { SAMPLE_SUMMARY } from '../../modules/local/report/sample_summary.nf'
include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'
include { FILTER_SUMMARY } from '../../modules/local/report/filter_summary.nf'

workflow GENERATE_REPORT {
    take:
    vcf_pre
    tbi_pre
    vcf_post
    tbi_post
    vcf_select
    tbi_select
    cv_file
    snpio_pre
    snpio_post
    clumpp_pre
    clumpp_post
    inds
    pops
    metrics
    top_loci

    main:
    ch_versions = Channel.empty()
    ch_mqc_files = Channel.empty()

    //CV plot
    PLOT_CV( cv_file )
    ch_versions = ch_versions.mix( PLOT_CV.out.versions )
    ch_mqc_files = ch_mqc_files.mix( PLOT_CV.out.cv_html )

    //Get individual lists from vcfs
    BCFTOOLS_QUERY_PRE( vcf_pre, tbi_pre )
    BCFTOOLS_QUERY_POST( vcf_post, tbi_post )
    ch_versions = ch_versions.mix( BCFTOOLS_QUERY_PRE.out.versions )

    //SNPio summary
    SAMPLE_SUMMARY(
        BCFTOOLS_QUERY_PRE.out.samples,
        BCFTOOLS_QUERY_POST.out.samples,
        snpio_pre,
        snpio_post
    )
    ch_mqc_files = ch_mqc_files.mix( SAMPLE_SUMMARY.out.summary_txt )
    ch_versions = ch_versions.mix( SAMPLE_SUMMARY.out.versions )

    //Admixture comparison
    COMPARE_ADMIXTURE(
        clumpp_pre,
        clumpp_post,
        inds,
        pops
    )
    ch_mqc_files = ch_mqc_files
        | mix( COMPARE_ADMIXTURE.out.min_max_html )
        | mix( COMPARE_ADMIXTURE.out.entropy_html )
        | mix( COMPARE_ADMIXTURE.out.summary_json )
    ch_versions = ch_versions.mix( COMPARE_ADMIXTURE.out.versions )

    //Admixture barplots
    ch_clumpp_pre  = clumpp_pre.map  { meta, file -> tuple(meta + [id: 'pre'], file) }
    ch_clumpp_post = clumpp_post.map { meta, file -> tuple(meta + [id: 'post'], file) }
    ch_clumpp_all = ch_clumpp_pre.mix( ch_clumpp_post )
    ch_inds = inds.map { meta, path -> path }
    ch_pops = pops.map { meta, path -> path }
    ch_plot_admix = ch_clumpp_all
        .combine(ch_inds)
        .combine(ch_pops)
    PLOT_ADMIXTURE(
        ch_plot_admix
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE.out.admixture_html )
    ch_versions = ch_versions.mix( COMPARE_ADMIXTURE.out.versions )

    //SNPio plots
    FILTER_SUMMARY( snpio_pre )
    ch_mqc_files = ch_mqc_files.mix( FILTER_SUMMARY.out.sankey_html )

    //Locus ranking plots
    PLOT_METRICS( metrics, top_loci )
    ch_mqc_files = ch_mqc_files.mix( PLOT_METRICS.out.hist_html )
    ch_mqc_files = ch_mqc_files.mix( PLOT_METRICS.out.scatter_html )
    ch_versions = ch_versions.mix( PLOT_METRICS.out.versions )

    emit:
    mqc_files    = ch_mqc_files
    versions     = ch_versions
}
