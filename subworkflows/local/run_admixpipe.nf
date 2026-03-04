//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { ADMIXTUREPIPELINE } from '../../modules/local/admixpipe/admixturepipeline.nf'
include { CLUMPAK } from '../../modules/local/admixpipe/submitclumpak.nf'
include { CVSUM } from '../../modules/local/admixpipe/cvsum.nf'
include { DISTRUCT } from '../../modules/local/admixpipe/distructrerun.nf'
include { EVALADMIX } from '../../modules/local/admixpipe/evaladmix.nf'
include { BESTK } from '../../modules/local/bestK.nf'
include { EVANNO } from '../../modules/local/evanno.nf'

workflow RUN_ADMIXPIPE {
    take:
    vcf         // [ val(meta), *.vcf or *.vcf.gz ]
    ch_popmap   // [ val(meta), popmap file ]

    main:
    ch_versions = Channel.empty()

    // Branch input VCF by extension
    vcf
    | branch {
        vcfgz: it[1].name.endsWith('.vcf.gz')
        vcf:   it[1].name.endsWith('.vcf')
    }
    | set { ch_vcf_branch }

    // If input was vcf.gz, decompress
    TABIX_BGZIP( ch_vcf_branch.vcfgz )
    ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )

    // Combine uncompressed .vcf with bgzipped .vcf
    TABIX_BGZIP.out.output
        | mix( ch_vcf_branch.vcf )
        | set { ch_vcf }

    // Pass to ADMIXTURE pipeline
    ADMIXTUREPIPELINE(
        ch_vcf,
        ch_popmap
    )
    ch_versions = ch_versions.mix( ADMIXTUREPIPELINE.out.versions )

    // Run CLUMPAK
    CLUMPAK(
        ADMIXTUREPIPELINE.out.results,
        ADMIXTUREPIPELINE.out.inds,
        ADMIXTUREPIPELINE.out.pops
    )
    ch_versions = ch_versions.mix( CLUMPAK.out.versions )

    // Run Distruct
    DISTRUCT(
        ADMIXTUREPIPELINE.out.pfiles,
        ADMIXTUREPIPELINE.out.qfiles,
        ADMIXTUREPIPELINE.out.pops,
        ADMIXTUREPIPELINE.out.inds,
        ADMIXTUREPIPELINE.out.logs,
        CLUMPAK.out.output
    )

    // Compute best K from crossval
    CVSUM(
        DISTRUCT.out.cv,
        DISTRUCT.out.loglik
    )
    ch_versions = ch_versions.mix( CVSUM.out.versions )

    // Run EvalAdmix
    EVALADMIX(
        ADMIXTUREPIPELINE.out.ped,
        ADMIXTUREPIPELINE.out.map,
        ADMIXTUREPIPELINE.out.pfiles,
        ADMIXTUREPIPELINE.out.qfiles,
        ADMIXTUREPIPELINE.out.qfiles_json,
        ch_popmap,
        CLUMPAK.out.output,
        DISTRUCT.out.major_clusters,
        DISTRUCT.out.cvruns_json,
        DISTRUCT.out.qfilepaths_json
    )
    ch_versions = ch_versions.mix( EVALADMIX.out.versions )

    // Evanno calculations (backup for bestK)
    EVANNO(
        CVSUM.out.ll_output
    )

    // Fetch results for the best K value
    BESTK(
        CVSUM.out.cv_output,
        DISTRUCT.out.best_results
    )
    ch_versions = ch_versions.mix( BESTK.out.versions )

    emit:
    best_results = DISTRUCT.out.best_results
    bestK_file   = BESTK.out.bestK_file
    bestK_clumpp = BESTK.out.bestK_clumpp
    evanno       = EVANNO.out.metrics
    inds         = ADMIXTUREPIPELINE.out.inds
    pops         = ADMIXTUREPIPELINE.out.pops
    cv_file      = CVSUM.out.cv_output
    qfilepaths   = DISTRUCT.out.qfilepaths_json
    corres       = EVALADMIX.out.corres
    fam          = EVALADMIX.out.fam
    versions     = ch_versions
}
