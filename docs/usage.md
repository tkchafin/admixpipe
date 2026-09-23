# aCaMEL/admixpipe: Usage

## Introduction

aCaMEL/admixpipe estimates population structure from a multi-sample VCF using ADMIXTURE. This page covers the input files, every pipeline parameter, how the best K is chosen, how to add map layers, and how to run and tune the pipeline. For a description of the results, see [output.md](output.md).

## Inputs

### Genotypes (`--input`)

One multi-sample VCF, uncompressed (`.vcf`) or bgzip-compressed (`.vcf.gz`). The pipeline creates the tabix index itself.

- The VCF should contain biallelic SNPs. Multi-allelic sites and indels should be removed beforehand, e.g. with `bcftools view -m2 -M2 -v snps`.
- Contig names and positions are used for physical thinning (`--thin_dist`). Data mapped to a pseudo-reference (e.g. concatenated RAD loci) works, provided each locus has its own contig or loci are separated by at least `--thin_dist` bp.
- The file's base name (everything before `.vcf`) is used as the prefix for output files.

BCF input is not currently supported. Convert it with `bcftools view -Oz -o out.vcf.gz in.bcf`.

### Population map (`--popmap`)

A tab-delimited text file with **no header** and two columns: sample ID, then population (or sampling site) ID.

```text
86NCCFE01	NCCFE
86NCCFE02	NCCFE
86NCCOT01	NCCOT
```

- Sample IDs must match the VCF header exactly.
- Every sample in the VCF should be listed in the popmap. With SNPio's `force_popmap` behaviour, samples missing from either file are dropped, so check the per-sample table in the report to confirm who was retained.
- Population IDs are used to group and order individuals in barplots, for per-population missing-data filters and summaries, for F<sub>ST</sub>, and to join samples to coordinates in `--site_coords`.

### Site coordinates (`--site_coords`, optional)

A tab-delimited text file with **no header** and three columns: population/site ID (matching the popmap), latitude, longitude, in decimal degrees (WGS84).

```text
NCCFE	38.37291666666618	-96.49375000000065
NCCOT	38.38541666666615	-96.54791666666732
```

The first column may instead hold **individual** IDs, one row per sample, to map samples collected at different points. The pipeline matches rows against population IDs and against sample IDs, and uses whichever matches more rows.

When this file is supplied, the report adds interactive maps of mean ancestry per site (pie charts, sized by sample count), for the best K and for every K. Samples with no matching coordinates are left off the maps.

### Map layers (`--geo_data_config`, `--geo_data_dir`, optional)

Extra vector layers, such as rivers, watersheds or range boundaries, can be drawn beneath the pie charts. Two parameters are needed:

- `--geo_data_dir`: a directory holding the layer files. Use any format GeoPandas can read, e.g. GeoPackage, GeoJSON or a shapefile with its sidecar files.
- `--geo_data_config`: a JSON list describing each layer.

```json
[
  {
    "path": "my_layers/streams.gpkg",
    "z_order": 1,
    "style": { "color": "#0066ff", "weight": 1, "opacity": 1.0 }
  },
  {
    "path": "my_layers/range.geojson",
    "z_order": 0,
    "style": { "color": "#555555", "weight": 2, "fillOpacity": 0.1 }
  }
]
```

- `path` is resolved relative to the **parent** of `--geo_data_dir`, so it starts with the directory's own name. With `--geo_data_dir /data/my_layers`, the first layer above is read from `/data/my_layers/streams.gpkg`.
- `z_order` sets the drawing order (lower values are drawn first, underneath).
- `style` is passed to Leaflet as a [path style](https://leafletjs.com/reference.html#path-option).

The bundled example is [`assets/test_geo_data.json`](../assets/test_geo_data.json) with [`assets/test_geo_data/`](../assets/test_geo_data/). Map layers require `--site_coords`.

## Parameters

Every parameter can also be listed with `nextflow run aCaMEL/admixpipe --help`. Every value used in a run, defaults included, is recorded in the _Workflow Summary_ section of the report and in `pipeline_info/params_<timestamp>.json`.

### Input/output

| Parameter           | Default   | Description                                                                 |
| ------------------- | --------- | --------------------------------------------------------------------------- |
| `--input`           | required  | Input VCF (`.vcf` or `.vcf.gz`).                                            |
| `--popmap`          | required  | Population map (see above).                                                 |
| `--outdir`          | `results` | Output directory. Use an absolute path on cloud storage.                    |
| `--site_coords`     | none      | Site coordinates. Enables the map sections of the report.                   |
| `--geo_data_config` | none      | JSON describing extra map layers. Requires `--geo_data_dir`.                |
| `--geo_data_dir`    | none      | Directory containing the files referenced by `--geo_data_config`.           |
| `--multiqc_title`   | none      | Title shown at the top of the report and used in its file name.             |
| `--email`           | none      | Address to send a completion summary to.                                    |

### Filtering

Filters are applied by SNPio in this order: per-population SNP missingness (`--pop_cov`), monomorphic sites, physical thinning (`--thin_dist`), overall SNP missingness (`--snp_cov`), individual missingness (`--ind_cov`), and minor allele frequency (`--min_maf`). The _Summary of Filtering Steps_ Sankey diagram in the report shows how many loci each step removed.

| Parameter     | Default | Description                                                                                     |
| ------------- | ------- | ----------------------------------------------------------------------------------------------- |
| `--ind_cov`   | `0.9`   | Maximum proportion of missing genotypes allowed per individual.                                 |
| `--snp_cov`   | `0.9`   | Maximum proportion of missing genotypes allowed per SNP.                                        |
| `--pop_cov`   | `0.9`   | Maximum proportion of missing genotypes allowed per SNP within any population.                  |
| `--min_maf`   | `0.05`  | Minimum minor allele frequency.                                                                 |
| `--thin_dist` | `100`   | SNPs within this many bp of another SNP are removed, to reduce linkage between retained SNPs.   |

ADMIXTURE assumes unlinked loci, so keep thinning on for data with many SNPs per locus. For RADseq-style data, a `--thin_dist` longer than the read length keeps roughly one SNP per locus.

### ADMIXTURE and choice of K

| Parameter        | Default | Description                                                            |
| ---------------- | ------- | ---------------------------------------------------------------------- |
| `--maxk`         | `10`    | Largest K tested. ADMIXTURE is run for every K from 1 to `--maxk`.     |
| `--num_reps`     | `10`    | Independent ADMIXTURE replicates per K (different random seeds).       |
| `--num_cv`       | `10`    | Number of cross-validation folds.                                      |
| `--bestk_method` | `cv`    | How the best K is chosen: `cv`, `evanno`, `lnl`, `l1` or `l2`.         |

The total number of ADMIXTURE runs is `maxk × num_reps`, so runtime grows with both values.

#### Choosing the best K

All of the results below are reported whichever method you choose. The only effect of `--bestk_method` is which K is highlighted and used for the _Best K_ sections of the report. K = 1 is never selected.

| Method   | Chooses the K that…                                                                                        |
| -------- | ---------------------------------------------------------------------------------------------------------- |
| `cv`     | has the lowest mean ADMIXTURE cross-validation error across replicates.                                    |
| `evanno` | has the highest Evanno ΔK = mean(\|L″(K)\|) / sd(L(K)) (Evanno _et al._ 2005).                              |
| `lnl`    | sits at the "elbow" of the mean log-likelihood curve L(K).                                                 |
| `l1`     | sits at the elbow of the first-order rate of change L′(K).                                                 |
| `l2`     | sits at the elbow of \|L″(K)\|.                                                                            |

The elbow is the K furthest from the straight line joining the first and last points of the curve. Ties go to the smaller K.

No single criterion is reliable in every case. ΔK, for example, cannot select K = 1 and tends to favour K = 2 under hierarchical structure. Treat the best K as a starting point. Look at the barplots for all K, the evalAdmix residuals, and the PCA before interpreting results.

### Hidden and advanced options

The standard nf-core options for resource limits (`--max_cpus`, `--max_memory`, `--max_time`), institutional configs, notifications and MultiQC customisation are also available. Show them with `--help --validationShowHiddenParams`.

`--kriging` (interpolated ancestry surfaces) is present but not currently supported. Setting it stops the pipeline with an error.

## Running the pipeline

A typical command:

```bash
nextflow run aCaMEL/admixpipe \
    -profile docker \
    --input genotypes.vcf.gz \
    --popmap popmap.tsv \
    --outdir results
```

The pipeline will create the following files in your working directory:

```bash
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in the specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other Nextflow hidden files, e.g. history of pipeline runs and old logs.
```

### Parameter files

To reuse the same settings, put them in a YAML or JSON file and pass it with `-params-file`:

```bash
nextflow run aCaMEL/admixpipe -profile docker -params-file params.yaml
```

```yaml
input: "genotypes.vcf.gz"
popmap: "popmap.tsv"
site_coords: "site_coords.tsv"
outdir: "results"
maxk: 12
num_reps: 20
bestk_method: "evanno"
```

> [!WARNING]
> Do not use `-c <file>` to set parameters. Custom config files given with `-c` should only be used for resource requests, infrastructure settings, or module arguments (`ext.args`, see below).

### Updating the pipeline

Nextflow caches pipeline code on first run and keeps using the cached version. To get the latest version:

```bash
nextflow pull aCaMEL/admixpipe
```

### Reproducibility

Specify a release with `-r` (e.g. `-r 1.0.0`) so that the same code and containers are used every time. The version, commit ID, full command line and every parameter value are recorded in the report's _Workflow Summary_ section. A draft methods paragraph with references is in its _Methods Description_ section. Keep `pipeline_info/params_<timestamp>.json` and any custom config files with your results.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double hyphen).

### `-profile`

Selects configuration presets. Several can be combined, e.g. `-profile test,docker`. Later profiles override earlier ones.

A container engine is required. Several steps use purpose-built containers for AdmixPipe, SNPio and plotting, so the `conda` and `mamba` profiles do **not** work with this pipeline.

- `docker`: use [Docker](https://docker.com/).
- `singularity` / `apptainer`: use [Singularity](https://sylabs.io/docs/) or [Apptainer](https://apptainer.org/). Recommended on HPC systems.
- `podman`, `shifter`, `charliecloud`: use those engines.
- `arm`: add this alongside `docker` on Apple Silicon or other ARM machines. It runs the x86-64 images under emulation, e.g. `-profile docker,arm`.
- `test`: runs a small bundled dataset (with maps and a map layer) and needs no other parameters.
- `test_full`: the same dataset with the default ADMIXTURE settings (K = 1–10, 10 replicates, 10-fold CV).

The pipeline also loads institutional profiles from [nf-core/configs](https://github.com/nf-core/configs), so `-profile <your_institution>` may already work on your cluster.

### `-resume`

Restarts a run, reusing cached results for any step whose inputs have not changed. This is useful after adjusting a report-only setting, or after a failure late in the run. You can resume a specific run with `-resume <run-name>`. `nextflow log` lists previous run names.

### `-c`

Loads an extra config file, e.g. for resource requests or tool arguments (below).

## Custom configuration

### Resource requests

Each step has default CPU, memory and time requests (see [`conf/base.config`](../conf/base.config)). Steps that fail with a resource-related exit code are retried automatically with larger requests. The upper limits are set by `--max_cpus`, `--max_memory` and `--max_time`.

The ADMIXTURE step does most of the work. It uses up to 12 CPUs and has a 48 h time limit, doubled on retry. To change these, add something like this to a config file and pass it with `-c`:

```groovy
process {
    withName: 'ADMIXTUREPIPELINE' {
        cpus   = 24
        memory = 32.GB
        time   = 96.h
    }
}
```

### Tool arguments

Extra arguments can be passed to individual steps with `ext.args` in a config file. For example, to change the colour palette and basemap of the ancestry maps:

```groovy
process {
    withName: 'PLOT_ADMIXTURE_SPATIAL' {
        ext.args = '--palette Set2 --basemap CartoDB.Positron'
    }
}
```

Some steps already set `ext.args` in [`conf/modules.config`](../conf/modules.config). If you override one, include the existing arguments you want to keep.

Setting `ext.args` for `SNPIO_FILTER` to `--permutations 1000` gives permutation p-values for pairwise F<sub>ST</sub>. The default is 0, i.e. no permutations. This can add substantially to runtime.

### Running on HPC

Use `-profile singularity` or `-profile apptainer`, plus an executor config for your scheduler, e.g.:

```groovy
process.executor = 'slurm'
process.queue    = 'comp'
```

Or use your institution's profile from nf-core/configs if one exists.

## Troubleshooting

- **Samples missing from results**: check that the sample IDs in `--popmap` exactly match the VCF header (`bcftools query -l file.vcf.gz`). Also check the _Missing Data Per-sample_ section of the report, since samples above `--ind_cov` are removed.
- **Very few SNPs after filtering**: check the _Summary of Filtering Steps_ Sankey diagram to see which filter removed most loci. Relax that filter, e.g. raise `--pop_cov` for populations with few samples, or lower `--thin_dist`.
- **Maps are missing**: maps are only made when `--site_coords` is given, and only for populations present in that file.
- **Maps have no background tiles**: basemap tiles load from the internet when the report is opened, so view it with a network connection.
- **Apple Silicon**: add the `arm` profile, e.g. `-profile docker,arm`.

## Running in the background

Nextflow must keep running until the pipeline finishes. Use `nextflow run ... -bg` to detach it from the terminal, or run it inside `screen` or `tmux`. On HPC you can also submit Nextflow itself as a (long, low-resource) job.

## Nextflow memory requirements

To stop the Nextflow Java process itself from using too much memory, add this to your environment, e.g. in `~/.bashrc`:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
