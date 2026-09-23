# aCaMEL/admixpipe: Output

## Introduction

This page describes the files the pipeline writes and the sections of its HTML report. Paths are relative to the directory given with `--outdir`. `<prefix>` is the base name of the input VCF, and `<prefix>_filtered` is the prefix used for results computed on the filtered data.

Most people only need the report:

```text
<outdir>/report/multiqc_report.html
```

The report is self-contained and can be shared as a single file. Map basemap tiles load from the internet when it is opened.

## Pipeline overview

- [Report](#report): the interactive report and the figures it contains
- [Filtering (SNPio)](#filtering-snpio): the filtered VCF and SNPio summaries
- [ADMIXTURE runs](#admixture-runs): raw Q/P matrices and logs for every replicate
- [CLUMPAK](#clumpak): replicate alignment and clustering modes
- [distruct](#distruct): aligned major-mode results and barplots per K
- [Cross-validation and log-likelihood summaries](#cross-validation-and-log-likelihood-summaries)
- [Evanno metrics](#evanno-metrics)
- [Best K selection](#best-k-selection)
- [evalAdmix](#evaladmix): model-fit residuals
- [Pipeline information](#pipeline-information): run metadata for reproducibility

## Report

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `multiqc_report.html`: the interactive report. Open it in a web browser.
  - `multiqc_data/`: data behind the report, as parsed by MultiQC.
  - `multiqc_plots/`: static exports of the report's standard plots.
  - `*_mqc.html`, `*_mqc.json`: the individual report sections as stand-alone files.
  - `admixture_spatial.tsv`: mean ancestry proportion per site at the best K, with coordinates and sample counts (only with `--site_coords`).
  - `*.samples.txt`: sample names in the input and filtered VCFs.
  - `geo_data_files/`: staged copies of map layers (only with `--geo_data_config`).

</details>

The report is built with [MultiQC](http://multiqc.info) and has these sections, in order.

### Dataset Composition (SNPio)

- **Summary of Filtering Steps**: a Sankey diagram of how many loci each filter removed and how many were kept. Use it to see which filter is limiting your data.
- **Missing Data Per-sample**: missingness for every sample before and after filtering, and whether it was retained.
- **Missing Data Per-population**: the same, summarised per population.
- **Pairwise Fst (θst)**: a heatmap of Weir & Cockerham F<sub>ST</sub> between populations on the filtered data.
- **PCA Scatter Plot**: the first principal components of the filtered genotypes, coloured by population. This gives a model-free view of structure to compare against ADMIXTURE.

### Best K

These sections show the K chosen by `--bestk_method`, which is marked in every plot.

- **Cross-validation Error Across K**: mean ± spread of ADMIXTURE CV error per K, across replicates.
- **Evanno Method (ΔK) Across K**: L(K), L′(K), \|L″(K)\| and ΔK.
- **Barplot**: individual ancestry proportions at the best K, from the CLUMPAK-aligned major mode, grouped by population.
- **Ancestry Proportions Map** (with `--site_coords`): an interactive map with one pie chart per site showing mean ancestry, plus any extra map layers.
- **EvalAdmix**: the residual correlation matrix for the best K. Values near zero mean the model fits. Positive correlations within a population mean that population's structure is not captured by the model (e.g. K too low, or a history ADMIXTURE cannot represent). Negative values between populations can mean spurious shared ancestry.

### All Results

The same plots for every K from 2 to `--maxk`, with a selector to switch between K values:

- **ADMIXTURE Assignment Probability Barplot**: barplots for all K.
- **Ancestry Proportions Map** (with `--site_coords`): maps for all K.
- **EvalAdmix**: residual correlation matrices for all K.

### Workflow Summary

Records everything needed to reproduce the run:

- **Command line**: the exact `nextflow run ...` command.
- **Run information**: pipeline version, git revision and commit ID, run name, session ID, start time, Nextflow version, profiles, config files, container engine and directories.
- **Parameters**: every pipeline parameter with the value used and its default. Values changed from the default are in bold. Standard nf-core options are in collapsed sections, and any extra parameters not in the pipeline schema are listed separately.

### Software Versions

The version of each tool used, per step.

### Methods Description

A draft methods paragraph for this run, listing the tools used with in-text citations and a reference list, plus the command line.

## Filtering (SNPio)

<details markdown="1">
<summary>Output files</summary>

- `snpio/snpio_filter/`
  - `<prefix>.filter.vcf.gz`, `<prefix>.filter.vcf.gz.tbi`: the filtered VCF used for all downstream steps.
  - `<prefix>_output/`: SNPio results for the unfiltered data (missingness reports and plots, logs).
  - `filtered_output/`: SNPio results for the filtered data, including population-genetic summary statistics, F<sub>ST</sub> and PCA.
  - `multiqc/multiqc_report.html`: SNPio's own stand-alone report, with more detail than the main report (allele frequency spectra, heterozygosity, locus-level missingness).
  - `multiqc/multiqc_report_data/`: tables behind the SNPio report.

</details>

[SNPio](https://github.com/btmartin721/SNPio) filters loci and samples as described in [usage.md](usage.md#filtering), and computes the data summaries shown in the _Dataset Composition_ section of the report.

## ADMIXTURE runs

<details markdown="1">
<summary>Output files</summary>

- `admixturepipeline/`
  - `<prefix>_filtered.<K>_<rep>.Q`: ancestry proportions (individuals × K) for each replicate.
  - `<prefix>_filtered.<K>_<rep>.P`: allele frequencies (SNPs × K) for each replicate.
  - `<prefix>_filtered.<K>_<rep>.stdout`: ADMIXTURE log for each run, including CV error and log-likelihood.
  - `<prefix>_filtered.ped`, `<prefix>_filtered.map`: the PLINK files ADMIXTURE was run on.
  - `<prefix>_filtered_inds.txt`, `<prefix>_filtered_pops.txt`: individual and population order of the rows in the Q files.
  - `<prefix>_filtered.qfiles.json`: an index of the Q files.
  - `results.zip`: all runs, packaged for CLUMPAK.

</details>

[ADMIXTURE](https://dalexander.github.io/admixture/) is run through [AdmixPipe](https://github.com/stevemussmann/admixturePipeline) for K = 1 to `--maxk`, with `--num_reps` replicates per K. Rows of the Q files follow the order in `*_inds.txt`.

## CLUMPAK

<details markdown="1">
<summary>Output files</summary>

- `clumpak/clumpakOutput/`
  - `K=<K>/`: CLUMPAK results for each K, including the major and any minor clustering modes and their CLUMPP-aligned Q matrices.
  - `K=<K>.MajorCluster.png`: barplot of the major mode for each K.

</details>

[CLUMPAK](https://clumpak.tau.ac.il/) groups the replicate runs at each K into clustering modes and aligns cluster labels across replicates. If a K has minor modes, the replicates did not all converge on the same solution. This is worth checking before interpreting that K.

## distruct

<details markdown="1">
<summary>Output files</summary>

- `distruct/`
  - `<prefix>_filtered/best_results/`
    - `ClumppIndFile.output.<K>`: aligned individual ancestry proportions for the major mode at each K. The report's barplots and maps are built from these.
    - `ClumppPopFile.<K>`: the same, averaged per population.
    - `K<K>.ps`: distruct barplots (PostScript).
  - `*.pdf`: the distruct barplots as PDFs.
  - `MajorClusterRuns.txt`: which replicates belong to the major mode at each K.
  - `cv_file.MajClust.txt`, `loglikelihood_file.MajClust.txt`: CV error and log-likelihood of the major-mode replicates.
  - `cvRuns.json`, `qfilePaths.json`: indexes used by later steps.

</details>

## Cross-validation and log-likelihood summaries

<details markdown="1">
<summary>Output files</summary>

- `cvsum/`
  - `cv_output.txt`: CV error summary per K (mean and spread across major-mode replicates).
  - `ll_output.txt`: log-likelihood summary per K.
  - `cv_file.MajClust.png`, `loglikelihood_file.MajClust.png`: static plots of the above.

</details>

## Evanno metrics

<details markdown="1">
<summary>Output files</summary>

- `evanno/`
  - `evanno_metrics.tsv`: per K, the mean and SD of the log-likelihood, L′(K), L″(K) and ΔK (Evanno _et al._ 2005).

</details>

## Best K selection

<details markdown="1">
<summary>Output files</summary>

- `bestk/`
  - `bestK.txt`: the K selected by `--bestk_method`.
  - `best_clumpp_indfile.out`: aligned individual ancestry proportions at that K.

</details>

## evalAdmix

<details markdown="1">
<summary>Output files</summary>

- `evaladmix/`
  - `*.corres`: residual correlation matrices from [evalAdmix](https://github.com/GenisGE/evalAdmix), for each K.
  - `*.png`: static plots of the residual correlations for the major (and any minor) modes and for each replicate.
  - `*.fam`: sample order of the correlation matrices.

</details>

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `params_<timestamp>.json`: every parameter value used, in machine-readable form. Pass it back with `-params-file` to re-run with the same settings.
  - `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt`, `pipeline_dag_<timestamp>.html`: Nextflow's resource-usage report, timeline, per-task trace and workflow graph.
  - `nf_core_pipeline_software_mqc_versions.yml`: software versions of every tool used.
  - `pipeline_report.html`, `pipeline_report.txt`: the completion summary, written only when `--email` / `--email_on_fail` is set.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) produces its own reports for troubleshooting and for estimating resource needs. Together with the report's _Workflow Summary_ section, these files record how a run was performed.
