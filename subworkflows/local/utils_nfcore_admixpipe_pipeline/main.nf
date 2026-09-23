//
// Subworkflow with functionality specific to the pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip/main'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir            // string: The output directory where the results will be saved
    input             // string: Path to input VCF or VCF.gz file
    popmap            // string: path to popmap file
    site_coords
    geo_data_config
    geo_data_dir

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input input.vcf[.gz] --popmap popmap.tsv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromPath(input)
        .map { file ->
            def meta = [id: file.simpleName]
            return [meta, file]
        }
        .branch {
            vcf: it[1].name.endsWith('.vcf')
            vcfgz: it[1].name.endsWith('.vcf.gz')
        }
        .set { ch_input }

    // Process VCF inputs
    TABIX_BGZIP ( ch_input.vcf )
    ch_tabix_vcf_input = ch_input.vcfgz
        | mix (TABIX_BGZIP.out.output )
    TABIX_TABIX( ch_tabix_vcf_input )

    //
    // Create channel for the popmap
    //
    Channel
        .fromPath(popmap)
        .map { file ->
            def meta = [id: file.simpleName]
            return [meta, file]
        }
        .set{ ch_popmap }

    //
    // Channel for geo_data_config (optional)
    //
    if ( params.geo_data_config ) {
        Channel
            .fromPath( params.geo_data_config )
            .map { file ->
                def meta = [ id: file.simpleName ]
                return [ meta, file ]
            }
            .set { ch_geo_data_config }
    }
    else {
        Channel
            .empty()
            .set { ch_geo_data_config }
    }


    //
    // Channel for a *pre‑staged* geodata directory (optional)
    //
    if ( params.geo_data_dir ) {
        Channel
            .fromPath( params.geo_data_dir )      // accepts dir or wildcard
            .map { dir ->
                def meta = [ id: file(dir).getBaseName() ]
                return [ meta, dir ]
            }
            .set { ch_geo_data_dir }
    }
    else {
        Channel.empty().set { ch_geo_data_dir }
    }

    //
    // Channel for site_coords (optional)
    //
    if ( params.site_coords ) {
        Channel
            .fromPath( params.site_coords )
            .map { file ->
                def meta = [ id: file.simpleName ]
                return [ meta, file ]
            }
            .set { ch_site_coords }
    }
    else {
        Channel
            .empty()
            .set { ch_site_coords }
    }

    // Collect versions
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf       = ch_tabix_vcf_input
    tbi       = TABIX_TABIX.out.tbi
    popmap    = ch_popmap
    site_coords = ch_site_coords
    geo_data    = ch_geo_data_config
    geo_data_dir = ch_geo_data_dir
    versions  = ch_versions
}


/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // -------------------------
    // maxk
    // -------------------------
    if (!(params.maxk instanceof Integer)) {
        try {
            params.maxk = params.maxk as Integer
        } catch (Exception e) {
            error("Invalid value for --maxk: '${params.maxk}'. It must be an integer.")
        }
    }

    if (params.maxk <= 1) {
        error("Invalid value for --maxk: '${params.maxk}'. It must be > 1.")
    }

    // -------------------------
    // Coverage parameters
    // -------------------------
    ['ind_cov','snp_cov','pop_cov','min_maf'].each { p ->
        if (!(params[p] instanceof Number)) {
            try {
                params[p] = params[p] as Double
            } catch (Exception e) {
                error("Invalid value for --${p}: '${params[p]}'. It must be numeric.")
            }
        }

        if (params[p] < 0 || params[p] > 1) {
            error("Invalid value for --${p}: '${params[p]}'. Must be between 0 and 1.")
        }
    }

    // -------------------------
    // thin_dist
    // -------------------------
    if (!(params.thin_dist instanceof Integer)) {
        try {
            params.thin_dist = params.thin_dist as Integer
        } catch (Exception e) {
            error("Invalid value for --thin_dist: '${params.thin_dist}'. It must be an integer.")
        }
    }

    // -------------------------
    // num_cv
    // -------------------------
    if (!(params.num_cv instanceof Integer)) {
        try {
            params.num_cv = params.num_cv as Integer
        } catch (Exception e) {
            error("Invalid value for --num_cv: '${params.num_cv}'. It must be an integer.")
        }
    }

    // -------------------------
    // num_reps
    // -------------------------
    if (!(params.num_reps instanceof Integer)) {
        try {
            params.num_reps = params.num_reps as Integer
        } catch (Exception e) {
            error("Invalid value for --num_reps: '${params.num_reps}'. It must be an integer.")
        }
    }

    // -------------------------
    // kriging (not supported yet)
    // -------------------------
    if (params.kriging) {
        error("The --kriging option is not currently supported. Please set --kriging false.")
    }

    // -------------------------
    // geodata layers
    // -------------------------
    if (params.geo_data_config && !params.geo_data_dir) {
        error("--geo_data_config requires --geo_data_dir (the directory holding the layer files).")
    }

    if (params.geo_data_config && !params.site_coords) {
        log.warn "--geo_data_config is ignored because --site_coords was not provided."
    }

    // -------------------------
    // bestk_method
    // -------------------------
    def allowedBestK = ["cv","lnl","l1","l2","evanno"]

    if (params.bestk_method == null) {
        error("Parameter --bestk_method cannot be null.")
    }

    params.bestk_method = params.bestk_method.toString().toLowerCase()

    if (!allowedBestK.contains(params.bestk_method)) {
        error("Invalid value for --bestk_method: '${params.bestk_method}'. Allowed values: ${allowedBestK.join(', ')}")
    }
}

//
// Full parameter summary for MultiQC
//
// Unlike paramsSummaryMultiqc (which only lists values that differ from the
// schema defaults), this lists every schema parameter alongside its default,
// any parameters passed that are not in the schema, and the exact command
// line and run configuration.
//
def fullParamsSummaryMultiqc(schema_filename) {
    def schema = new groovy.json.JsonSlurper().parse(file("${workflow.projectDir}/${schema_filename}"))

    def esc = { v ->
        v == null ? '' : v.toString()
            .replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            .replaceAll(/\r?\n/, ' ')
    }
    def fmt = { v ->
        (v == null || v.toString() == '') ? '<span style="color:#999999;">N/A</span>' : "<samp>${esc(v)}</samp>"
    }
    def table = { List header, List rows ->
        def out = ['<table class="table table-condensed table-hover" style="width:auto;">']
        out << '<thead><tr>' + header.collect { "<th>${it}</th>" }.join('') + '</tr></thead><tbody>'
        rows.each { row -> out << '<tr>' + row.collect { "<td>${it}</td>" }.join('') + '</tr>' }
        out << '</tbody></table>'
        return out
    }

    def lines = []

    // Run information
    lines << '<h4>Command line</h4>'
    lines << "<pre><code>${esc(workflow.commandLine)}</code></pre>"
    lines << '<h4>Run information</h4>'
    def run_info = [
        'Pipeline version' : workflow.manifest.version,
        'Revision'         : workflow.revision,
        'Commit ID'        : workflow.commitId,
        'Run name'         : workflow.runName,
        'Session ID'       : workflow.sessionId,
        'Started'          : workflow.start,
        'Nextflow version' : workflow.nextflow.version,
        'Profile'          : workflow.profile,
        'Config files'     : workflow.configFiles?.join(', '),
        'Container engine' : workflow.containerEngine,
        'Launch directory' : workflow.launchDir,
        'Work directory'   : workflow.workDir,
        'Project directory': workflow.projectDir,
        'User'             : workflow.userName
    ]
    lines.addAll(table(['Field', 'Value'], run_info.collect { k, v -> ["<b>${k}</b>", fmt(v)] }))

    // Parameters, grouped as in the schema
    lines << '<h4>Parameters</h4>'
    lines << '<p>Values that differ from the pipeline default are shown in <b>bold</b>.</p>'
    def seen = [] as Set
    schema.definitions.each { group_id, group ->
        def props = group['properties'] ?: [:]
        if (!props) return
        def rows = props.collect { name, spec ->
            seen << name
            def value   = params.containsKey(name) ? params[name] : null
            def changed = value?.toString() != spec['default']?.toString()
            def label   = changed ? "<b>${name}</b>" : name
            [label, changed ? "<b>${fmt(value)}</b>" : fmt(value), fmt(spec['default'])]
        }
        // Pipeline-specific groups are expanded; nf-core boilerplate groups are collapsed
        def boilerplate = group_id in ['institutional_config_options', 'max_job_request_options', 'generic_options']
        lines << (boilerplate ? "<details><summary><b>${esc(group.title)}</b></summary>" : "<p style=\"font-size:110%\"><b>${esc(group.title)}</b></p>")
        lines.addAll(table(['Parameter', 'Value', 'Default'], rows))
        if (boilerplate) lines << '</details>'
    }

    // Anything passed that the schema does not know about (e.g. from -params-file)
    def extra = params.keySet().findAll { !(it in seen) && !it.contains('-') }.sort()
    if (extra) {
        lines << '<p style="font-size:110%"><b>Other parameters (not in schema)</b></p>'
        lines.addAll(table(['Parameter', 'Value'], extra.collect { [it, fmt(params[it])] }))
    }

    String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
    yaml_file_text        += "description: ' - full command line, run configuration and parameter values for this run.'\n"
    yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text        += "section_href: '${workflow.manifest.homePage}'\n"
    yaml_file_text        += "plot_type: 'html'\n"
    yaml_file_text        += "data: |\n"
    yaml_file_text        += lines.collect { "    ${it}" }.join('\n') + '\n'

    return yaml_file_text
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Input genotypes were compressed and indexed with tabix (Li 2011), and sample lists were extracted with bcftools (Danecek et al. 2021).",
            "SNPs and individuals were filtered, and missingness, F<sub>ST</sub> and PCA summaries were computed, with SNPio (Martin et al. 2026).",
            "Ancestry proportions were estimated with ADMIXTURE (Alexander et al. 2009) via AdmixPipe (Mussmann et al. 2020, 2023), using VCFtools (Danecek et al. 2011) and PLINK (Chang et al. 2015) for file conversion.",
            "Replicate runs were aligned with CLUMPAK (Kopelman et al. 2015) and CLUMPP (Jakobsson & Rosenberg 2007), and plotted with distruct (Rosenberg 2004).",
            "Model fit was assessed with evalAdmix (Garcia-Erill & Albrechtsen 2020), and the best K was chosen using ${params.bestk_method == 'evanno' ? 'the Evanno delta K method (Evanno et al. 2005)' : params.bestk_method == 'cv' ? 'ADMIXTURE cross-validation error' : 'the ' + params.bestk_method + ' criterion (Evanno et al. 2005)'}.",
            "Results were summarised with MultiQC (Ewels et al. 2016)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19(9), 1655–1664. doi: <a href='https://doi.org/10.1101/gr.094052.109'>10.1101/gr.094052.109</a></li>",
            "<li>Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4, 7. doi: <a href='https://doi.org/10.1186/s13742-015-0047-8'>10.1186/s13742-015-0047-8</a></li>",
            "<li>Danecek, P., Auton, A., Abecasis, G., et al. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156–2158. doi: <a href='https://doi.org/10.1093/bioinformatics/btr330'>10.1093/bioinformatics/btr330</a></li>",
            "<li>Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: <a href='https://doi.org/10.1093/gigascience/giab008'>10.1093/gigascience/giab008</a></li>",
            "<li>Evanno, G., Regnaut, S., & Goudet, J. (2005). Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Molecular Ecology, 14(8), 2611–2620. doi: <a href='https://doi.org/10.1111/j.1365-294X.2005.02553.x'>10.1111/j.1365-294X.2005.02553.x</a></li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href='https://doi.org/10.1093/bioinformatics/btw354'>10.1093/bioinformatics/btw354</a></li>",
            "<li>Garcia-Erill, G., & Albrechtsen, A. (2020). Evaluation of model fit of inferred admixture proportions. Molecular Ecology Resources, 20(4), 936–949. doi: <a href='https://doi.org/10.1111/1755-0998.13171'>10.1111/1755-0998.13171</a></li>",
            "<li>Jakobsson, M., & Rosenberg, N. A. (2007). CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics, 23(14), 1801–1806. doi: <a href='https://doi.org/10.1093/bioinformatics/btm233'>10.1093/bioinformatics/btm233</a></li>",
            "<li>Kopelman, N. M., Mayzel, J., Jakobsson, M., Rosenberg, N. A., & Mayrose, I. (2015). Clumpak: a program for identifying clustering modes and packaging population structure inferences across K. Molecular Ecology Resources, 15(5), 1179–1191. doi: <a href='https://doi.org/10.1111/1755-0998.12387'>10.1111/1755-0998.12387</a></li>",
            "<li>Li, H. (2011). Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5), 718–719. doi: <a href='https://doi.org/10.1093/bioinformatics/btq671'>10.1093/bioinformatics/btq671</a></li>",
            "<li>Martin, B. T., Monaco, D. R., Sharabi, N., Mussmann, S. M., & Chafin, T. K. (2026). SNPio: a Python interface for population genomic data processing. BMC Bioinformatics. doi: <a href='https://doi.org/10.1186/s12859-026-06546-5'>10.1186/s12859-026-06546-5</a></li>",
            "<li>Mussmann, S. M., Douglas, M. R., Chafin, T. K., & Douglas, M. E. (2020). AdmixPipe: population analyses in Admixture for non-model organisms. BMC Bioinformatics, 21, 337. doi: <a href='https://doi.org/10.1186/s12859-020-03701-4'>10.1186/s12859-020-03701-4</a></li>",
            "<li>Mussmann, S. M., Douglas, M. R., Chafin, T. K., & Douglas, M. E. (2023). AdmixPipe v3: facilitating population structure delimitation from SNP data. Bioinformatics Advances, 3(1), vbad168. doi: <a href='https://doi.org/10.1093/bioadv/vbad168'>10.1093/bioadv/vbad168</a></li>",
            "<li>Rosenberg, N. A. (2004). distruct: a program for the graphical display of population structure. Molecular Ecology Notes, 4(1), 137–138. doi: <a href='https://doi.org/10.1046/j.1471-8286.2003.00566.x'>10.1046/j.1471-8286.2003.00566.x</a></li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
