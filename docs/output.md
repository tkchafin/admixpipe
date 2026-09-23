# aCaMEL/admixpipe: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. `<prefix>` is the base name of the input VCF.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [SNPio](#snpio) - Filtered VCF, missing data, F<sub>ST</sub> and PCA
- [ADMIXTURE](#admixture) - Ancestry estimates for every K and replicate
- [CLUMPAK](#clumpak) - Replicate alignment and clustering modes
- [distruct](#distruct) - Aligned results and barplots for every K
- [Choice of K](#choice-of-k) - Cross-validation, log-likelihood and Evanno summaries, and the best K
- [evalAdmix](#evaladmix) - Model fit for every K
- [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### SNPio

<details markdown="1">
<summary>Output files</summary>

- `snpio/snpio_filter/`
  - `<prefix>.filter.vcf.gz`: filtered VCF used for all later steps, with its `.tbi` index.
  - `<prefix>_output/`, `filtered_output/`: SNPio results before and after filtering.
  - `multiqc/multiqc_report.html`: SNPio's own report, with more detail than the main report.

</details>

[SNPio](https://github.com/btmartin721/SNPio) filters SNPs and samples by missing data, minor allele frequency and physical distance. It also summarises missing data per sample and per population, and computes pairwise F<sub>ST</sub> and a PCA.

### ADMIXTURE

<details markdown="1">
<summary>Output files</summary>

- `admixturepipeline/`
  - `<prefix>_filtered.<K>_<rep>.Q`: ancestry proportions for each K and replicate.
  - `<prefix>_filtered.<K>_<rep>.P`: allele frequencies for each K and replicate.
  - `<prefix>_filtered.<K>_<rep>.stdout`: ADMIXTURE log, including CV error and log-likelihood.
  - `<prefix>_filtered_inds.txt`, `<prefix>_filtered_pops.txt`: sample and population order of the Q files.

</details>

[ADMIXTURE](https://dalexander.github.io/admixture/) is run by [AdmixPipe](https://github.com/stevemussmann/admixturePipeline) for K = 1 to `--maxk`, with `--num_reps` replicates each.

### CLUMPAK

<details markdown="1">
<summary>Output files</summary>

- `clumpak/clumpakOutput/`
  - `K=<K>/`: clustering modes and aligned results for each K.
  - `K=<K>.MajorCluster.png`: barplot of the major mode for each K.

</details>

[CLUMPAK](https://clumpak.tau.ac.il/) aligns replicate runs and groups them into modes. If a K has minor modes, its replicates did not all converge on the same solution.

### distruct

<details markdown="1">
<summary>Output files</summary>

- `distruct/`
  - `<prefix>_filtered/best_results/ClumppIndFile.output.<K>`: aligned ancestry proportions for the major mode at each K.
  - `*.pdf`: barplots for each K.
  - `MajorClusterRuns.txt`: replicates in the major mode at each K.

</details>

[distruct](https://rosenberglab.stanford.edu/distruct.html) plots the aligned major-mode results for each K.

### Choice of K

<details markdown="1">
<summary>Output files</summary>

- `cvsum/`
  - `cv_output.txt`, `ll_output.txt`: cross-validation error and log-likelihood for each K.
- `evanno/`
  - `evanno_metrics.tsv`: L(K), L′(K), L″(K) and Evanno ΔK for each K.
- `bestk/`
  - `bestK.txt`: the K chosen by `--bestk_method`.

</details>

### evalAdmix

<details markdown="1">
<summary>Output files</summary>

- `evaladmix/`
  - `*.corres`: residual correlation matrix for each K and replicate.
  - `*.png`: plots of the residual correlations.

</details>

[evalAdmix](https://github.com/GenisGE/evalAdmix) tests how well the ADMIXTURE model fits the data. Residual correlations near zero mean a good fit.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.
  - `*_mqc.html`, `*_mqc.json`: individual report sections.
  - `admixture_spatial.tsv`: mean ancestry per site at the best K (with `--site_coords`).

</details>

[MultiQC](http://multiqc.info) collects the results into one interactive report, with these sections:

- **Dataset Composition**: filtering summary, missing data per sample and population, pairwise F<sub>ST</sub> and PCA.
- **Best K**: CV error and Evanno plots, the barplot, map and evalAdmix residuals at the best K.
- **All Results**: barplots, maps and evalAdmix residuals for every K.
- **Workflow Summary**: the full command line, run information and every parameter value, with defaults.
- **Software Versions** and **Methods Description**: tool versions, and a methods paragraph with references.

Maps are only included when `--site_coords` is given. Map backgrounds are loaded from the internet when the report is opened.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and `pipeline_dag_<timestamp>.html`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `nf_core_pipeline_software_mqc_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameters are used when running the pipeline.
  - Parameters used by the pipeline run: `params_<timestamp>.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
