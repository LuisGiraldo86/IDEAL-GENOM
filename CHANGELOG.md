# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.2] - 2026-09-02

### Fixed
- `ideal_genom.__version__` always reported `"0.0.0"`. `__init__.py` looked
  up `importlib.metadata.version("ideal-genom-qc")`, the package's PyPI name
  before it was renamed to `ideal-genom`; that lookup always raised
  `PackageNotFoundError`. Now looks up `"ideal-genom"`.

### Added
- `.github/workflows/release.yml` — pushing a `vX.Y.Z` tag now builds,
  verifies the tag matches `pyproject.toml`'s version, publishes to PyPI via
  Trusted Publishing (no stored token), and creates a GitHub Release with the
  built sdist/wheel attached. Replaces manual `poetry build`/`poetry publish`.

## [1.3.1] - 2026-09-02

### Fixed
- **`ref_annotation` validated before the post-imputation pipeline runs, not
  after.** `AnnotateVCF`/`ProcessVCF.execute_process_vcf_pipeline()` passed
  `--annotations <file> --columns ID` to `bcftools annotate` without checking
  that the reference file was actually BGZF-compressed and tabix/CSI-indexed.
  A plain-gzip or unindexed reference file made bcftools silently fall back to
  its tab-annotation mode, failing with an unrelated `The -c CHROM option not
  given` error — only surfaced at the very last step of the pipeline, after
  unzip/filter/normalize/index had already run. Added
  `validate_bgzip_indexed_vcf()` in `ideal_genom/core/utils.py`, called from
  `AnnotateVCF.__init__` and at the top of
  `ProcessVCF.execute_process_vcf_pipeline()`, so a bad `ref_annotation` now
  fails immediately with an actionable message.

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
