include { PLOT_CV } from '../../modules/local/report/plot_cv.nf'
include { PLOT_EVANNO } from '../../modules/local/report/plot_evanno.nf'
include { SAMPLE_SUMMARY } from '../../modules/local/report/sample_summary.nf'
include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'
include { PLOT_ADMIXTURE_KRIGING } from '../../modules/local/report/plot_admixture_kriging.nf'
include { PLOT_ADMIXTURE_KRIGING_MULTIK } from '../../modules/local/report/plot_admixture_kriging_multik.nf'
include { PLOT_ADMIXTURE_MULTIK } from '../../modules/local/report/plot_admixture_all.nf'
include { PLOT_ADMIXTURE_SPATIAL } from '../../modules/local/report/plot_admixture_spatial.nf'
include { PLOT_ADMIXTURE_SPATIAL_MULTIK } from '../../modules/local/report/plot_admixture_spatial_multik.nf'
include { FILTER_SUMMARY } from '../../modules/local/report/filter_summary.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_PRE } from '../../modules/local/bcftools_query.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_POST } from '../../modules/local/bcftools_query.nf'
include { POP_SUMMARY } from '../../modules/local/report/pop_summary.nf'
include { PLOT_PAIRWISE_FST } from '../../modules/local/report/pairwise_fst.nf'
include { PLOT_PCA } from '../../modules/local/report/plot_pca.nf'
include { PLOT_EVALADMIX } from '../../modules/local/report/plot_evaladmix.nf'
include { STAGE_GEODATA_LAYERS } from '../../modules/local/report/stage_geodata.nf'

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
    snpio_report_data
    clumpp
    inds
    pops
    qfilepaths
    corres
    fam
    site_coords
    geo_data
    geo_data_dir

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


    //Admixture barplots
    PLOT_ADMIXTURE(
        clumpp,
        inds,
        pops
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE.out.versions )

    //Admixture barplots -- all
    PLOT_ADMIXTURE_MULTIK(
        best_results,
        inds,
        pops,
        bestk_file
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_MULTIK.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE_MULTIK.out.versions )

    //Admixture maps
    if (params.site_coords){
        if (params.geo_data_config){

            STAGE_GEODATA_LAYERS( geo_data, geo_data_dir )

            PLOT_ADMIXTURE_SPATIAL(
                clumpp,
                inds,
                pops,
                site_coords,
                STAGE_GEODATA_LAYERS.out.geo_data_dir
            )
            ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_SPATIAL.out.plot_html )
            ch_versions = ch_versions.mix( PLOT_ADMIXTURE_SPATIAL.out.versions )

            if (params.kriging){
                PLOT_ADMIXTURE_KRIGING(
                    clumpp,
                    inds,
                    pops,
                    site_coords,
                    STAGE_GEODATA_LAYERS.out.geo_data_dir
                )
                ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_KRIGING.out.html_discrete )
                    .mix( PLOT_ADMIXTURE_KRIGING.out.html_simpson )
                ch_versions = ch_versions.mix( PLOT_ADMIXTURE_KRIGING.out.versions )

                PLOT_ADMIXTURE_KRIGING_MULTIK(
                    best_results,
                    inds,
                    pops,
                    site_coords,
                    STAGE_GEODATA_LAYERS.out.geo_data_dir
                )
            }

            //ADMIXTURE maps (all K)
            PLOT_ADMIXTURE_SPATIAL_MULTIK(
                best_results,
                inds,
                pops,
                site_coords,
                STAGE_GEODATA_LAYERS.out.geo_data_dir
            )
            ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_SPATIAL_MULTIK.out.plot_html )
            ch_versions = ch_versions.mix( PLOT_ADMIXTURE_SPATIAL_MULTIK.out.versions )

        }else{
            PLOT_ADMIXTURE_SPATIAL(
                clumpp,
                inds,
                pops,
                site_coords,
                tuple( [], [] )
            )
            ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_SPATIAL.out.plot_html )
            ch_versions = ch_versions.mix( PLOT_ADMIXTURE_SPATIAL.out.versions )

            if (params.kriging){
                PLOT_ADMIXTURE_KRIGING(
                    clumpp,
                    inds,
                    pops,
                    site_coords,
                    tuple( [], [] )
                )
                ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_KRIGING.out.html_discrete )
                    .mix( PLOT_ADMIXTURE_KRIGING.out.html_simpson )
                ch_versions = ch_versions.mix( PLOT_ADMIXTURE_KRIGING.out.versions )

                PLOT_ADMIXTURE_KRIGING_MULTIK(
                    best_results,
                    inds,
                    pops,
                    site_coords,
                    tuple( [], [] )
                )
            }

            PLOT_ADMIXTURE_SPATIAL_MULTIK(
                best_results,
                inds,
                pops,
                site_coords,
                tuple( [], [] )
            )
            ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE_SPATIAL_MULTIK.out.plot_html )
            ch_versions = ch_versions.mix( PLOT_ADMIXTURE_SPATIAL_MULTIK.out.versions )

        }
    }

    //EvalAdmix (best and all K)
    PLOT_EVALADMIX( qfilepaths, fam, corres, bestk_file)
    ch_mqc_files = ch_mqc_files.mix( PLOT_EVALADMIX.out.allk_html )
    ch_mqc_files = ch_mqc_files.mix( PLOT_EVALADMIX.out.bestk_html )

    //Get individual lists from vcfs
    BCFTOOLS_QUERY_PRE( vcf_pre, tbi_pre )
    BCFTOOLS_QUERY_POST( vcf_post, tbi_post )
    ch_versions = ch_versions.mix( BCFTOOLS_QUERY_PRE.out.versions )

    //SNPio sample missingness
    SAMPLE_SUMMARY(
        snpio_report_data,
        BCFTOOLS_QUERY_PRE.out.samples,
        BCFTOOLS_QUERY_POST.out.samples
    )
    ch_mqc_files = ch_mqc_files.mix( SAMPLE_SUMMARY.out.summary_txt )
    ch_versions = ch_versions.mix( SAMPLE_SUMMARY.out.versions )

    //SNPio sample missingness
    POP_SUMMARY(
        snpio_report_data
    )
    ch_mqc_files = ch_mqc_files.mix( POP_SUMMARY.out.summary_txt )
    ch_versions = ch_versions.mix( POP_SUMMARY.out.versions )

    //SNPio sample missingness
    PLOT_PAIRWISE_FST(
        snpio_report_data
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_PAIRWISE_FST.out.plot_html )
    ch_versions = ch_versions.mix( PLOT_PAIRWISE_FST.out.versions )

    //SNPio sankey
    FILTER_SUMMARY( snpio_pre )
    ch_mqc_files = ch_mqc_files.mix( FILTER_SUMMARY.out.sankey_html )

    //SNPio PCA
    PLOT_PCA(
        snpio_report_data
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_PCA.out.plot_html )
    ch_versions = ch_versions.mix( PLOT_PCA.out.versions )

    emit:
    mqc_files    = ch_mqc_files
    versions     = ch_versions
}
