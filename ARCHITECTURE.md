# Architecture: control & data flow

This document traces how a run actually executes, end-to-end, through the two entry points into IDEAL-GENOM. It complements `CLAUDE.md` (which maps modules to responsibilities) by showing how those pieces call each other at runtime.

## Two entry points, one set of step classes

```
                    ┌─────────────────────────┐
   YAML config ───▶ │   CLI / PipelineExecutor │   "production" path: orchestrated, chained via config
                    └───────────┬─────────────┘
                                │
                                ▼
                  ┌──────────────────────────┐
   direct calls ─▶│   Step classes            │◀── notebooks/*.ipynb  "exploratory" path: manual, one
   (no orchestr.) │   (SampleQC, AncestryQC,   │                       notebook cell per call, manual chaining
                  │   VariantQC, GWAS_GLM, …) │
                  └───────────┬──────────────┘
                              │
                              ▼
                  ┌──────────────────────────┐
                  │  core/executor.py         │
                  │  run_plink / run_plink2   │──▶ PLINK1.9 / PLINK2 / GCTA / bcftools subprocess
                  └──────────────────────────┘
```

Both paths bottom out in the same step classes and the same PLINK/GCTA/bcftools subprocess calls — the orchestrator only adds config parsing, step-chaining, and automatic report/cleanup invocation on top.

## Path A: CLI (`ideal-genom run --config pipeline.yaml`)

```
cli.py cmd_run
  └─ config.load_config(path)        # parses YAML, validates schema, raises ConfigurationError on failure
  └─ PipelineExecutor(config).execute()
       │
       ├─ Pass 1 — instantiate every step (including disabled ones)
       │     for step_config in pipeline.steps:
       │         _instantiate_step(step_config)   # so later steps can reference even a disabled step's attrs
       │         self.steps[name] = instance
       │
       └─ Pass 2 — for each *enabled* step, in YAML order:
             _resolve_params(init_params / execute_params)
               • "${base_output_dir}"          → self.base_output_dir
               • "${steps.sample_qc.clean_dir}" → getattr-chain into self.steps['sample_qc']
             _convert_paths_to_path_objects()  # any *_path / *_file string → pathlib.Path
             instance = module.Class(**init_params)         # e.g. SampleQC(...)
             instance.execute_<step_name>_pipeline(execute_params)   # runs the whole stage, see "Inside a stage" below
             self.steps[step_name] = instance
             _generate_reports(step_name, instance)   # if settings.reports.generate_reports
                 → instantiates the matching *Report class, feeds it instance's output-file attributes
             _perform_cleanup(step_name, instance)    # if not settings.files.keep_intermediate
                 → instantiates the matching *CleanUp class (sample_qc / variant_qc only)
```

Step-chaining is what makes `ancestry_qc`'s `input_path: "${steps.sample_qc.clean_dir}"` work in the YAML templates: by the time `ancestry_qc` is reached, `self.steps['sample_qc']` already holds the fully-executed `SampleQC` instance, and its `clean_dir` attribute is read directly off that live object — there's no intermediate serialization.

## Path B: direct Python API (notebooks)

The notebooks in `notebooks/` do by hand what `PipelineExecutor` automates:

```python
step = SampleQC(input_path=..., input_name=..., output_path=..., output_name=..., high_ld_regions_file=...)
step.execute_sample_qc_pipeline(sample_params)        # same method PipelineExecutor would call

report = SampleQCReport(output_path=step.plots_dir)    # caller builds *Report manually
report.report_sample_qc(call_rate_smiss=step.call_rate_miss, ...)   # using attributes the step set on itself

cleanup = SampleQCCleanUp(output_path=step.results_dir, input_path=step.input_path)
cleanup.clean_all()
```

Chaining across stages is manual too — `02-ancestry_qc.ipynb` points its `input_path` directly at `outputData/sample_qc_results/clean_files/`, the literal path `01-sample_qc.ipynb` wrote to, rather than reading a live attribute off a `SampleQC` instance.

## Inside a stage: what `execute_<name>_pipeline()` actually does

Every stage's all-in-one method (`SampleQC.execute_sample_qc_pipeline`, `VariantQC.execute_variant_qc_pipeline`, `AncestryQC.execute_ancestry_qc_pipeline`, …) is a fixed sequence of `execute_*()` steps, and each of those follows the same inner loop:

```
execute_<check>()
  └─ core.utils.get_optimal_threads() / get_available_memory()   # size --threads / --memory from host
  └─ core.executor.run_plink(...) / run_plink2(...)              # subprocess call, writes files under self.results_dir
       (e.g. --missing, --check-sex, --het, --genome, --king-cutoff, --pca, ...)
  └─ self.<some_attr> = path to the PLINK output file just written
       (e.g. self.call_rate_miss, self.sexcheck_miss, self.maf_greater_het, self.kinship_miss, self.eigenvectors)
```

After all per-check steps have run, `get_fail_<samples|variants>()` reads those output files back into pandas (chunked, for large files), applies the configured thresholds, concatenates the failures, writes `fail_*.txt`/`fail_summary.txt`/`variant_qc_summary.tsv` to `self.fails_dir` / `self.results_dir`, and `execute_drop_<samples|variants|ancestry_outliers>()` does one final PLINK `--remove`/`--exclude` call against `self.input_path`, writing the cleaned `.bed/.bim/.fam` to `self.clean_dir`. That `clean_dir` output is what the next stage (whether chained by `PipelineExecutor` or copy-pasted by hand into a notebook) treats as its input.

## Cross-cutting: auto-fetched reference data

Any constructor that needs a reference file it wasn't given (`high_ld_regions_file`, 1000 Genomes `reference_files`) checks for it eagerly and fetches on miss, *before* any PLINK call happens:

```
__init__(..., high_ld_regions_file=path_that_may_not_exist):
    if not high_ld_regions_file.is_file():
        FetcherLDRegions(build=self.build).get_ld_regions()   # downloads into ideal_genom/data/ld_regions_files/
        high_ld_regions_file = fetcher.ld_regions              # now points at the cached download
```

This recurs verbatim (same shape, different fetcher) in `SampleQC`, `AncestryQC`, `ReferenceGenomicMerger`, `PCAReduction` (population/projection.py), and the GWAS `Preparatory` class — `core/get_references.py`'s `FetcherLDRegions`/`Fetcher1000Genome` are the shared implementation. Downloaded files are cached in `ideal_genom/data/`, so only the first run per build/region-set actually hits the network.
