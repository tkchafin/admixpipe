# aCaMEL/admixpipe: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of aCaMEL/admixpipe.

### `Added`

- VCF filtering and data summaries (missingness, pairwise F<sub>ST</sub>, PCA, filtering Sankey) with SNPio.
- ADMIXTURE across K = 1…`--maxk` with replicates and cross-validation via AdmixPipe, with CLUMPAK/CLUMPP replicate alignment and distruct plots.
- Best-K selection by cross-validation, Evanno ΔK, or the elbow of L(K), L′(K) or |L″(K)| (`--bestk_method`).
- Model-fit assessment with evalAdmix for all K.
- Interactive MultiQC report with barplots, CV/Evanno plots and evalAdmix heatmaps for the best K and all K.
- Optional interactive ancestry maps (`--site_coords`) with user-supplied vector layers (`--geo_data_config`, `--geo_data_dir`).
- Full run provenance in the report: command line, run metadata, and every parameter value with its default.
- Draft methods text with tool citations in the report.
