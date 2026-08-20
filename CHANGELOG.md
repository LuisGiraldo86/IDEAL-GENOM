# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-08-20

### Added
- **Brisbane plots** — `brisbane_draw()` and `brisbane_process_data()` in
  `ideal_genom/visualizations/manhattan_type.py`. Bins variants into fixed-size
  genomic windows (`window_kb`, default 100) and plots per-window SNP density
  against genomic position, in the style of Yengo et al. (2022, Nature), with
  optional genome-wide mean/median density lines.
- `get_yengo_height_independent_signals()` in `ideal_genom/core/get_examples.py`
  — downloads the 12,111 conditionally independent COJO signals from the GIANT
  height GWAS (Yengo et al. 2022, Supplementary Table 5) as example data.
  Coordinates are hg19/GRCh37.
- `viz_notebooks/brisbane.ipynb` — worked example for the Brisbane plot.

### Fixed
- **Post-imputation filtering no longer fails on Michigan Imputation Server
  output.** `FilterVariants` matched `unzipped-*.vcf.gz`, which picked up the
  per-chromosome `*.empiricalDose.vcf.gz` files alongside the intended
  `*.dose.vcf.gz` ones. Those files carry no `R2` INFO field, so
  `bcftools view -i 'R2>…'` aborted with exit code 255 and killed the pipeline.
  The pattern is now `unzipped-*.dose.vcf.gz`.
- **Chromosome ordering in Manhattan and Miami plots.** Cumulative chromosome
  offsets and row sorting used lexicographic order, placing chr10 before chr2
  and scattering X/Y/MT arbitrarily. Ordering is now genomic (1–22, X, Y, XY,
  MT/M, then unrecognized labels alphabetically), handling int, numeric-string
  and `chr`-prefixed labels. Existing Manhattan/Miami plots will change
  appearance where the input was not already in genomic order.
- Corrected an undefined-variable reference in `ProcessVCF.execute_concatenate()`
  that raised `NameError` instead of the intended `TypeError` when `output_name`
  was not a string.

### Changed
- `find_chromosomes_center()` and the Manhattan annotation helper rewritten as
  vectorized pandas operations instead of row-wise loops; behavior unchanged.
- Deduplicated redundant column-existence checks in `manhattan_draw()`.

[1.3.0]: https://github.com/LuisGiraldo86/IDEAL-GENOM/compare/v1.2.0...v1.3.0
