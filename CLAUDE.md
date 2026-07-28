# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

IDEAL-GENOM is a Python package (Poetry-managed, `ideal_genom/` package, distributed as `ideal-genom` on PyPI) for automated, reproducible genotype data analysis. It wraps PLINK1.9, PLINK2, GCTA and bcftools as subprocesses and provides three pipelines: sample/ancestry/variant **QC**, post-imputation **VCF processing**, and **GWAS** (GLM/GLMM). It's usable as a Python API, via Jupyter notebooks (`notebooks/`), or via a YAML-driven CLI.

## Commands

Install for development:
```bash
poetry install
```

Run the CLI (entry point `ideal-genom = ideal_genom.__main__:main`, defined in `pyproject.toml`):
```bash
ideal-genom run --config my_pipeline.yaml [--dry-run]
ideal-genom validate --config my_pipeline.yaml
ideal-genom template --output my_pipeline.yaml [--force]
ideal-genom --version
```
Example YAML configs to start from: `yaml_configs/qc_pipeline_config_template.yaml`, `gwas_config_template.yaml`, `vcf_config_template.yaml`.

Build docs (Sphinx, config at `docs/source/conf.py`):
```bash
cd docs && make html
```

Build/run the Docker image (PLINK1.9, PLINK2, GCTA, bcftools preinstalled):
```bash
docker build -t ideal-genom .
docker run -it ideal-genom
```

**There is no test suite** — `tests/` contains only an empty `__init__.py`, and no pytest/ruff/mypy config exists anywhere in the repo. Verify changes by running the relevant pipeline class directly, via the CLI against a YAML config, or by exercising the matching notebook in `notebooks/`.

## Architecture

### Pipeline-stage class pattern

Every analysis stage (sample QC, variant QC, ancestry QC, GWAS GLM/GLMM, VCF processing, dimensionality reduction) follows the same three-way split:

- **`Execute*` / main class** (e.g. `SampleQC`, `VariantQC`, `AncestryQC`, `GWAS_GLM`) — does the actual work by shelling out to PLINK/GCTA/bcftools (via `core/executor.py`'s `run_plink`/`run_plink2`), with one `execute_<step>()` method per logical step plus a `execute_<name>_pipeline(params: dict)` convenience method that runs all steps in order. Constructors validate input paths/PLINK file triples (`.bed`/`.bim`/`.fam`) eagerly and raise `TypeError`/`FileNotFoundError`/`ValueError`.
- **`*Report` class** (e.g. `SampleQCReport`, `VariantQCReport`, `AncestryQCReport`) — separate class, takes file paths produced by the main class and only generates matplotlib/seaborn plots; no QC logic, no fail-list computation.
- **`*CleanUp` class** (e.g. `SampleQCCleanUp`, `VariantQCCleanUp`) — deletes intermediate PLINK byproducts. Not every module has one (ancestry QC and the dimensionality-reduction module don't; ancestry QC only has a private `_clean_merging_dir()` called automatically by its pipeline method).

When modifying one of these, check whether the "do the work" and "report on the work" responsibilities need to move together — they're intentionally decoupled, and the `*_pipeline()` convenience methods already call fail-aggregation and drop/cleanup steps internally, so don't re-implement that aggregation in a caller.

### Module map

- `ideal_genom/qc/` — `sample_qc.py`, `variant_qc.py`, `ancestry_qc.py`. Ancestry QC additionally contains `ReferenceGenomicMerger` (harmonizes study data with a 1000 Genomes reference: strand-ambiguous SNP filtering, chr/position/allele-flip fixes, merge) and `GenomicOutlierAnalyzer` (PCA + outlier detection by standardized distance from reference/study centroids — must exceed both thresholds to be flagged).
- `ideal_genom/population/projection.py` — population structure / dimensionality reduction: `PCAReduction` (LD pruning + PCA), `UMAPReduction`, `TSNEReduction` (pure `fit_transform()` on PCA eigenvectors), `Plot2D` (metadata-aware plotting, decoupled from the reduction step), and `DimensionalityReductionPipeline` orchestrating all of them. `execute_dimensionality_reduction_pipeline()` auto-detects whether `umap_params`/`tsne_params` contain list values vs. scalars and transparently switches between a single run and a full parameter-grid sweep (`execute_parameter_grid()`).
- `ideal_genom/population/fst_stats.py` — `FstSummary`: merges study data with a 1000 Genomes reference (via `ReferenceGenomicMerger`, reused from `qc/ancestry_qc.py`) and computes Fst statistics with PLINK, orchestrated by `execute_fst_pipeline(fst_params: dict)`.
- `ideal_genom/gwas/` — `preparatory.py` (`Preparatory`: LD pruning + PCA prep), `gen_linear_model.py` (`GWAS_GLM`: PLINK2 fixed-effects GLM), `gen_linear_mix_model.py` (`GWAS_GLMM`: GCTA GRM + mixed model), all with top-hits extraction and gene annotation.
- `ideal_genom/post_imputation/` — `vcf_process.py` (bcftools-based parallel pipeline: `UnzipVCF` → `FilterVariants` → `NormalizeVCF` → `ReferenceNormalizeVCF` → `IndexVCF` → `AnnotateVCF` (optional, only if an annotation reference is given) → concatenate, orchestrated by `ProcessVCF.execute_process_vcf_pipeline()`, with each step class built on a shared `ParallelTaskRunner` ThreadPoolExecutor base), `vcf_to_plink.py` (`GetPLINK`: VCF → PLINK1.9 binaries via `convert_vcf_to_plink()`).
- `ideal_genom/core/` — cross-cutting infrastructure:
  - `pipeline.py` — `PipelineExecutor`: the YAML-driven orchestrator. Dynamically imports `module`/`class` per step, resolves `${base_output_dir}` and `${steps.<name>.<attribute>}` references against already-instantiated step objects, converts any `*_path`/`*_file` string param to `Path`, and calls the step's `execute_<name>_pipeline()` method by naming convention. Disabled steps are still instantiated (not executed) so later steps can reference their attributes. Auto-generates reports/cleanup for `sample_qc`/`ancestry_qc`/`variant_qc` steps based on `settings.reports.generate_reports` / `settings.files.keep_intermediate`.
  - `cli.py` / `config.py` — `load_config()`/`validate_config()` enforce the YAML schema: `pipeline.{name,base_output_dir,steps}` required, each step needs `{name,module,class,init_params}` with `init_params.{input_path,input_name,output_path}` required (`execute_params`/`enabled` optional); raise `ConfigurationError` on violation. The top-level `settings` key is never validated and most of it is never read either — only `settings.files.keep_intermediate` and `settings.reports.{generate_reports,plot_format}` are actually consumed (by `PipelineExecutor`); `settings.logging`/`settings.resources` in the YAML templates are accepted but silently ignored — thread/memory sizing is always automatic via `core/utils.py`, not config-driven.
  - `executor.py` — `run_plink`/`run_plink2` subprocess wrappers used by every analysis class.
  - `utils.py` — `get_optimal_threads()`/`get_available_memory()` (auto-size PLINK `--threads`/`--memory` from host resources), `count_file_lines()`.
  - `get_references.py` — `FetcherLDRegions`, `Fetcher1000Genome`: auto-download high-LD-region files and 1000 Genomes reference panels into `ideal_genom/data/` when a class isn't given them explicitly (constructors check `is_file()`/`exists()` and fetch on miss — this pattern recurs in `SampleQC`, `AncestryQC`, `ReferenceGenomicMerger`, `PCAReduction`, GWAS `Preparatory`).
  - `environment.py` — `verify_genomic_environment()`/`quick_verify()` check PLINK/PLINK2/bcftools/GCTA are on `PATH` and report versions.
  - `get_examples.py` — downloads example GWAS summary-statistics datasets (BBJ, GCST cohorts) into `ideal_genom/data/sumstats/` for the `viz_notebooks/` tutorials.
- `ideal_genom/utilities/` — `annotations.py` (gene annotation via Ensembl/RefSeq GTF, `pyensembl`/`gtfparse`), `api_client.py` (`VEPEnsemblRestClient`, `GeneEnsemblRestClient`: rate-limited Ensembl REST clients with mirror fallback), `power_comp.py` (GWAS statistical power calculations).
- `ideal_genom/visualizations/` — `plots.py` (QQ plots, beta-beta plots, trumpet plots), `manhattan_type.py` (Manhattan/Miami plots), `zoom_heatmap.py`.

### Notebooks as canonical usage examples

`notebooks/01-sample_qc.ipynb`, `02-ancestry_qc.ipynb`, `03-variant_qc.ipynb`, `04-dimensionality_reduc.ipynb` are the maintained, runnable reference for the Python API (as opposed to the YAML/CLI path). They chain together: ancestry QC and dimensionality reduction read their input from the sample QC notebook's cleaned `PLINK` output (`outputData/sample_qc_results/clean_files/`). `viz_notebooks/` separately demonstrates the GWAS visualization functions using example summary statistics from `core/get_examples.py`. When changing a class's public method signatures, check these notebooks for breakage — they are not covered by any automated test.

### Genome build handling

Most classes take a `build: str` ('37' or '38') used only to pick which reference files (`FetcherLDRegions`/`Fetcher1000Genome`) to auto-fetch — it doesn't otherwise change behavior.
