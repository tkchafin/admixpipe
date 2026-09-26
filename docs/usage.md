# aCaMEL/admixpipe: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files. Run `nextflow run UARK-aCaMEL/admixpipe --help` to list them._

## Introduction

aCaMEL/admixpipe estimates population structure from SNP data with ADMIXTURE. The inputs are a VCF and a population map, plus optional site coordinates for maps.

## VCF input

A single multi-sample VCF of biallelic SNPs, as `.vcf` or `.vcf.gz`. The pipeline indexes it.

```bash
--input '[path to VCF file]'
```

## Population map input

A tab-delimited file with no header. The first column is the sample ID, which must match the VCF header. The second column is the population or site ID.

```bash
--popmap '[path to popmap file]'
```

```text title="popmap.tsv"
86NCCFE01	NCCFE
86NCCFE02	NCCFE
86NCCOT01	NCCOT
```

An [example popmap](../assets/test.popmap) is provided with the pipeline.

## Site coordinates input (optional)

A tab-delimited file with no header, with columns: population ID (or sample ID), latitude, longitude, in decimal degrees. When provided, the report includes interactive maps of ancestry proportions for each site.

```bash
--site_coords '[path to coordinates file]'
```

```text title="site_coords.tsv"
NCCFE	38.3729	-96.4938
NCCOT	38.3854	-96.5479
```

An [example coordinates file](../assets/site_coords.tsv) is provided with the pipeline.

### Map layers

Extra vector layers, such as rivers or range boundaries, can be drawn on the maps. Put the layer files in a directory, and describe them in a JSON file:

```bash
--geo_data_dir '[directory of layer files]' --geo_data_config '[layers JSON]'
```

```json title="layers.json"
[
  {
    "path": "my_layers/streams.gpkg",
    "z_order": 1,
    "style": { "color": "#0066ff", "weight": 1, "opacity": 1.0 }
  }
]
```

Each `path` starts with the name of the `--geo_data_dir` directory. `style` takes [Leaflet path options](https://leafletjs.com/reference.html#path-option). See [`assets/test_geo_data.json`](../assets/test_geo_data.json) for an example.

## Choosing the best K

`--bestk_method` sets which K is shown in the _Best K_ sections of the report. Results for every K are always reported.

| Method         | Best K is the one with                          |
| -------------- | ----------------------------------------------- |
| `cv` (default) | the lowest mean cross-validation error          |
| `evanno`       | the highest Evanno ΔK                           |
| `lnl`          | the elbow of the mean log-likelihood curve L(K) |
| `l1`           | the elbow of L′(K)                              |
| `l2`           | the elbow of \|L″(K)\|                          |

K = 1 is never selected. No single method is reliable in every case, so check the barplots, evalAdmix residuals and PCA before interpreting results.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run UARK-aCaMEL/admixpipe --input genotypes.vcf.gz --popmap popmap.tsv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run UARK-aCaMEL/admixpipe -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: 'genotypes.vcf.gz'
popmap: 'popmap.tsv'
outdir: './results/'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull UARK-aCaMEL/admixpipe
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [aCaMEL/admixpipe releases page](https://github.com/UARK-aCaMEL/admixpipe/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

The version number, the full command line and the value of every parameter are recorded in the _Workflow Summary_ section of the MultiQC report. The parameters are also saved to `pipeline_info/params_<timestamp>.json`, which can be passed back with `-params-file` to repeat a run.

> [!TIP]
> If you wish to share such a parameter file (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer) - see below.

> [!IMPORTANT]
> This pipeline requires a container engine. Conda is not supported, because several steps use purpose-built containers.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_full`
  - The same test data, run with the default ADMIXTURE settings
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `arm`
  - Use together with `docker` on Apple Silicon and other ARM machines
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code indicating it ran out of resources, it will automatically be resubmitted once with double the requests. If it fails again, the pipeline execution is stopped.

ADMIXTURE (`ADMIXTUREPIPELINE`) is the most demanding step. Its run time grows with `--maxk` × `--num_reps`.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
