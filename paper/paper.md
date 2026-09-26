---
title: 'aCaMEL/admixpipe: a reproducible Nextflow pipeline for ADMIXTURE-based population structure analysis'
tags:
  - Nextflow
  - nf-core
  - population genetics
  - population structure
  - ADMIXTURE
  - molecular ecology
  - conservation genomics
authors:
  - name: Tyler K. Chafin
    orcid: 0000-0001-8687-5905
    corresponding: true
    affiliation: 1
  - name: Steven M. Mussmann
    orcid: 0000-0002-5237-5088
    affiliation: 2
  - name: Marlis R. Douglas
    orcid: 0000-0001-6234-3939
    affiliation: 1
  - name: Michael E. Douglas
    orcid: 0000-0001-9670-7825
    affiliation: 1
affiliations:
  - name: Department of Biological Sciences, University of Arkansas, Fayetteville, Arkansas, United States
    index: 1
    ror: 05jbt9m15
  - name: Abernathy Fish Technology Center, U.S. Fish and Wildlife Service, Longview, Washington, United States
    index: 2
    ror: 04k7dar27
date: 24 September 2026
bibliography: paper.bib
---

# Summary

Most species are subdivided into partially isolated populations. Delineating these populations is foundational to evolutionary inference and to conservation, where genetic groups inform the designation of conservation and management units [@Funk2012; @Hohenlohe2021]. Model-based clustering has become the de facto standard for this task. It was popularised by STRUCTURE [@Pritchard2000; @Novembre2016] and extended to genome-scale data by ADMIXTURE [@Alexander2009]. These methods posit K ancestral populations with distinct allele frequencies and estimate each individual's ancestry as a mixture of them. A defensible analysis, however, is not a single model fit but a chain of dependent decisions. The data must be filtered, the model fitted repeatedly across K, replicate solutions reconciled, K chosen, and model adequacy assessed.

aCaMEL/admixpipe automates this chain end to end. From a VCF and a population map, it filters loci and individuals and runs replicated ADMIXTURE analyses. It then aligns replicates into clustering modes, applies alternative K-selection criteria, evaluates model fit and, optionally, projects ancestry onto maps. All results are consolidated into a single interactive report that also records the complete provenance of the analysis.

# Statement of need

AdmixPipe [@Mussmann2020] was developed because pipelines existed to summarise replicated STRUCTURE runs but none served ADMIXTURE, the faster maximum-likelihood alternative widely applied to reduced-representation data from non-model organisms. It coupled SNP filtering and replicated ADMIXTURE runs with CLUMPAK [@Kopelman2015]. AdmixPipe v3 [@Mussmann2023] added evalAdmix [@GarciaErill2020] tests within the multimodal context identified by CLUMPAK and containerised the software stack. Its authors characterised population structure inference as "a fundamental, but non-trivial task" compounded by "hierarchical population structure, diverse analytical methods, and complex software dependencies". AdmixPipe nonetheless executes as a monolithic process on a single host, leaving parallelisation, fault recovery and portability to the user.

Inference is also sensitive to choices made on either side of model fitting:

- minor allele frequency thresholds [@Linck2019] and missing-data filters [@HuangKnowles2016] alter the structure recovered;
- uneven sampling [@Puechmaille2016; @Wang2017] and hierarchical structure [@Kalinowski2011; @Janes2017] bias common K-selection criteria;
- barplots invite over-interpretation without formal assessment of fit [@Lawson2018].

These choices are rarely reported in full. A reanalysis of published STRUCTURE studies failed to recover the reported number of clusters in 30% of cases [@Gilbert2012].

These are not merely academic concerns. Genetic evidence increasingly informs the delineation of management units, listing decisions and stocking policy. Yet a persistent gap separates conservation genetics research from its implementation [@Taylor2017; @Kadykalo2020; @Klutsch2021]. A systematic review found definable management outcomes for only 49 of 115 applied studies [@Tkach2023]. When equivalent data are analysed under different, undocumented choices, results can appear contradictory, eroding their defensibility. Standardised workflows have been proposed as one bridge between science and practice [@Holderegger2019].

aCaMEL/admixpipe addresses these needs for population geneticists, molecular ecologists, and agency geneticists and their contractors. It belongs to a suite of standardised workflows for applied conservation genetics from the Arkansas Conservation and Molecular Ecology Lab (aCaMEL), alongside workflows for genotyping-panel design and hybrid classification.

# State of the field

Beyond STRUCTURE and ADMIXTURE [@Alexander2011], fastSTRUCTURE [@Raj2014] offers variational inference for large SNP datasets, and principal component analysis [@Patterson2006] provides a model-free complement. Post-processing remains fragmented across tools:

- STRUCTURE HARVESTER [@Earl2012] implements the Evanno method [@Evanno2005];
- CLUMPP [@Jakobsson2007] and CLUMPAK [@Kopelman2015] resolve label switching and multimodality among replicates;
- pong [@Behr2016] and pophelper [@Francis2017] provide visualisation;
- StructureSelector [@Li2018] aggregates alternative K estimators.

AdmixPipe integrates ADMIXTURE, CLUMPAK and evalAdmix but omits filtering diagnostics, spatial summaries and workflow management. scalepopgen [@Upadhyay2024], a general Nextflow toolkit for population genomics, fits ADMIXTURE once per K and plots cross-validation error and ancestry proportions. It does not align replicates, identify modes or assess fit.

Rather than reimplement this functionality, aCaMEL/admixpipe executes the published AdmixPipe container for ADMIXTURE, CLUMPAK, distruct [@Rosenberg2004] and evalAdmix, inheriting a validated core. Around it, the pipeline adds orchestration, filtering diagnostics, K selection, spatial summaries and provenance-rich reporting. Its authors include the developers of AdmixPipe, and filtering and summary statistics are provided by SNPio [@Martin2026].

# Software design

## Workflow framework

The pipeline is implemented in Nextflow DSL2 [@DiTommaso2017] on the nf-core template [@Ewels2020] (\autoref{fig:workflow}). Its cost scales with the number of K values multiplied by the number of replicates, and it consists largely of independent, long-running tasks. A workflow manager therefore provides parallel execution, checkpointed resumption and identical behaviour across workstations, clusters and cloud platforms [@Wratten2021]. Nextflow supports 18 schedulers or cloud services and seven container engines [@Langer2025].

nf-core contributes community-agreed best practice: schema-validated parameters, linting, continuous integration against bundled test data, and a library of reusable modules through which research communities can "adopt common standards progressively" [@Langer2025]. These conventions operationalise the FAIR principles for research software [@Wilkinson2016; @Barker2022], and they have measurable consequences. In an independent assessment, 51% of released nf-core pipeline revisions executed without failure, compared with 11% in the Snakemake Workflow Catalog. nf-core pipelines also remained reproducible for a median of 2.8 rather than 0.8 years [@Grayson2023].

![Overview of aCaMEL/admixpipe. The inset shows the steps executed within AdmixPipe.\label{fig:workflow}](figure1.png)

## Design decisions

- **Containers only.** Every step runs in a container [@Gruning2018]. Several dependencies, including the AdmixPipe stack, are distributed as purpose-built images rather than Conda packages. The pipeline therefore supports Docker, Singularity/Apptainer and Podman but not Conda, trading package-manager flexibility for an identical runtime environment.
- **Filtering in SNPio.** Delegating filtering to SNPio [@Martin2026] allows missingness per sample and per population to be reported before and after filtering, and locus attrition to be attributed to each filter. Pairwise F~ST~ and PCA are provided as model-free comparisons. Loci are physically thinned by default, consistent with ADMIXTURE's assumption of independence among loci.
- **K selection as a reported choice.** Five criteria are implemented: cross-validation error, Evanno ΔK, and the elbow of L(K), L′(K) or |L″(K)|. The chosen criterion determines only which K is emphasised. Barplots, evalAdmix residual correlations and maps are produced for every K, since no single criterion is reliable across demographic scenarios [@Janes2017; @Puechmaille2016].
- **Standardisation by construction.** Filters are applied in a fixed order under documented defaults. All results are consolidated into a single self-contained MultiQC [@Ewels2016] report, optionally with interactive maps of site-level ancestry over user-supplied vector layers. The report records the command line, pipeline revision, configuration profiles and every parameter value against its default, flagging deviations. It also generates a methods paragraph citing each tool used. It thereby automates the reporting recommended by @Gilbert2012, making analyses comparable across studies and auditable by those who commission them. Because the aCaMEL workflows share modules for filtering, ADMIXTURE execution and reporting, these conventions propagate across the suite.

# Research impact statement

aCaMEL/admixpipe builds on AdmixPipe, whose two descriptions have accrued 46 and 8 citations (Crossref, September 2026). The pipeline has been applied in agency-funded conservation assessment. It provided hybrid screening (K up to 20, with 20 replicates per K) and spatial population structure analyses for the endemic Beaded Darter (*Etheostoma clinton*), in a State Wildlife Grant report to the Arkansas Game and Fish Commission [@Bruckerhoff2026]. Its sibling workflow, aCaMEL/hybridclassification, shares its filtering, ADMIXTURE and reporting components. It underpinned hybrid classification in a genomic assessment of Smallmouth Bass for the same agency [@Douglas2026]. The release includes a bundled test dataset with sampling coordinates and a vector layer, continuous integration that executes the full workflow, and user and output documentation.

# AI usage disclosure

<!-- Authors: confirm this statement reflects all AI use across the project before submission. -->

Generative AI (Claude Opus 5.5, Anthropic, via Claude Code) assisted with this release. It was used to draft the user documentation and parameter schema, implement the report's parameter summary, identify and fix workflow bugs, verify references, and draft this paper. All AI-assisted changes were reviewed, edited and tested by the authors. The authors made all design decisions and take responsibility for the software and the paper.

# Acknowledgements

We thank Bradley T. Martin for SNPio and the nf-core community for the pipeline template. Applications of the pipeline were supported by the U.S. Fish and Wildlife Service State Wildlife Grants Program through the Arkansas Game and Fish Commission (AR-T-F22AF03392), and by the Arkansas Game and Fish Commission (SL4124). M.R.D. and M.E.D. acknowledge support from the Bruker Professorship in Life Sciences and the 21st Century Chair in Global Change Biology, respectively, at the University of Arkansas. The funders had no role in the design of the software or the preparation of this paper.

# References
