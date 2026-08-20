Changelog
=========

All notable changes to IDEAL-GENOM will be documented in this file.

Version 1.3.0 (Current)
-----------------------

**Released:** August 2026

**New Features:**

- **Brisbane plots** — ``brisbane_draw()`` and ``brisbane_process_data()`` in
  ``ideal_genom.visualizations.manhattan_type``. Bins variants into fixed-size genomic
  windows (``window_kb``, default 100) and plots per-window SNP density against genomic
  position, in the style of Yengo et al. (2022, *Nature*), with optional genome-wide
  mean and median density lines
- ``get_yengo_height_independent_signals()`` in ``ideal_genom.core.get_examples`` —
  downloads the 12,111 conditionally independent COJO signals from the GIANT height
  GWAS (Yengo et al. 2022, Supplementary Table 5) as example data. Coordinates are
  hg19/GRCh37
- ``viz_notebooks/brisbane.ipynb`` — worked example for the Brisbane plot

**Bug Fixes:**

- **Post-imputation filtering no longer fails on Michigan Imputation Server output.**
  ``FilterVariants`` matched ``unzipped-*.vcf.gz``, which picked up the per-chromosome
  ``*.empiricalDose.vcf.gz`` files alongside the intended ``*.dose.vcf.gz`` ones. Those
  files carry no ``R2`` INFO field, so ``bcftools view -i 'R2>…'`` aborted with exit
  code 255 and killed the pipeline. The pattern is now ``unzipped-*.dose.vcf.gz``
- **Chromosome ordering in Manhattan and Miami plots.** Cumulative chromosome offsets
  and row sorting used lexicographic order, placing chr10 before chr2 and scattering
  X/Y/MT arbitrarily. Ordering is now genomic (1–22, X, Y, XY, MT/M, then unrecognized
  labels alphabetically), handling int, numeric-string and ``chr``-prefixed labels.

  .. note::

     Existing Manhattan and Miami plots will change appearance where the input was not
     already in genomic order.

- Corrected an undefined-variable reference in ``ProcessVCF.execute_concatenate()`` that
  raised ``NameError`` instead of the intended ``TypeError`` when ``output_name`` was not
  a string

**Changed:**

- ``find_chromosomes_center()`` and the Manhattan annotation helper rewritten as
  vectorized pandas operations instead of row-wise loops; behavior unchanged
- Deduplicated redundant column-existence checks in ``manhattan_draw()``

**Documentation:**

- Documented the previously missing ``ideal_genom.core.environment``,
  ``ideal_genom.core.get_examples``, ``ideal_genom.utilities.api_client`` and
  ``ideal_genom.utilities.power_comp`` modules
- Added API reference entries for ``ReferenceGenomicMerger`` and
  ``GenomicOutlierAnalyzer``, the reference-harmonization and ancestry-outlier classes
  of the Ancestry QC module
- Added a Brisbane plot section to the visualization reference, based on
  ``viz_notebooks/brisbane.ipynb``
- Corrected the repository, issue-tracker and Read the Docs URLs throughout, and the
  ``pip install`` package name (``ideal-genom``, not ``ideal-genom-qc``)
- Corrected the notebook listing in ``examples.rst``, which named a nonexistent
  ``04-population.ipynb`` and omitted ``00-1KG_phenotype.ipynb`` and
  ``05-fst_statistic.ipynb``
- Removed a duplicated feature block from the documentation home page


Version 1.2.0
-------------

**Released:** July 2026

**Documentation:**

- Realigned the Sphinx docs with the actual codebase: fixed broken ``SampleQC``/``AncestryQC`` usage examples in ``api_overview.rst``, ``contributing.rst`` and ``examples.rst`` that referenced nonexistent constructor parameters and a nonexistent ``run_sample_qc()`` method
- Corrected ``configuration.rst`` to stop documenting ``settings.logging``/``settings.resources``/``settings.files`` keys that are accepted but silently ignored by the pipeline executor
- Rewrote ``visualization_modules.rst`` (manhattan_type, plots, zoom_heatmap) to use live ``automodule`` directives instead of hand-maintained function signatures that had drifted from the real code; usage examples rebuilt from the runnable ``viz_notebooks/`` notebooks
- Added the missing ``ideal_genom.population.fst_stats`` module to ``CLAUDE.md``'s module map
- Fixed stale ``0.2.0`` version references left over from the v0.2.0 YAML migration

**Bug Fixes:**

- Fixed ``docs/source/conf.py``'s ``autodoc_mock_imports``: it mocked ``matplotlib`` and other real project dependencies, which broke ``mpl_toolkits.mplot3d`` imports (via ``textalloc``) and pandas' own optional-``pyarrow`` version check during doc builds
- Fixed the ``ideal-genom --version`` fallback string (used only when package metadata lookup fails) still reading ``0.2.0``


Version 1.1.0
-------------

**Released:** January 2026

**New Features:**

- Modules for visualization of GWAS results
   - Manhattan plot
   - Miami plot
   - QQ plot
   - LD heatmap
   - Trumpet plot
   - Effect size comparison plot (beta-beta plot)
- Module to fetch example datasets for testing and tutorials

**Improvements:**

- Example notebooks for visualization modules


Version 1.0.0
-------------

**Released:** January 2026

**Major Changes:**

- Complete redesign of configuration system from JSON to YAML
- Unified pipeline framework with step-based execution
- New command-line interface with subcommands (run, validate, template)
- Package renamed from ideal-genom-qc to ideal-genom

**New Features:**

- **GWAS Pipeline**: Complete GWAS workflow with GLM and GLMM support
- **VCF Processing**: Post-imputation VCF processing and PLINK conversion
- **Population Analysis**: Enhanced population structure tools (PCA, UMAP, t-SNE, Fst)
- **Pipeline Executor**: Automated pipeline execution with dependency management
- **Configuration Validation**: Built-in YAML configuration validation
- **Variable Substitution**: Support for dynamic variable references in configurations

**Improvements:**

- Enhanced error handling and logging
- Improved memory management for large datasets
- Parallel processing support across all modules
- Better documentation with comprehensive examples
- Rich reporting and visualization capabilities
- Docker support with pre-configured environment

**API Changes:**

- YAML configuration replaces JSON configuration files
- New Python API with consistent class interfaces
- Module reorganization: ``ideal_genom.qc``, ``ideal_genom.gwas``, ``ideal_genom.post_imputation``
- Simplified class initialization patterns

**Bug Fixes:**

- Fixed kinship calculation memory issues
- Improved VCF parsing for large files
- Corrected reference genome download paths
- Fixed PCA projection edge cases

**Dependencies:**

- Python 3.11+ required (up from 3.8+)
- PLINK 2.0 now required alongside PLINK 1.9
- GCTA 1.95.0 or higher
- BCFtools 1.10+ for VCF processing

Version 0.1.0
-------------

**Released:** 2024

**Initial Release:**

- Basic QC pipeline (Sample QC, Ancestry QC, Variant QC)
- JSON-based configuration
- Command-line interface for basic operations
- Integration with PLINK 1.9
- 1000 Genomes reference panel support
- Basic UMAP visualization

Migration Guide (0.1.0 → 0.2.0)
--------------------------------

**Configuration Files:**

Old (JSON)::

    {
        "sample_qc": {
            "mind": 0.1,
            "maf": 0.01
        }
    }

New (YAML)::

    pipeline:
      steps:
        - name: "sample_qc"
          module: "ideal_genom.qc.sample_qc"
          class: "SampleQC"
          execute_params:
            mind: 0.1
            maf: 0.01

**Command-Line Interface:**

Old::

    python -m ideal_genom_qc \
        --path_params parameters.json \
        --file_folders paths.json \
        --steps steps.json

New::

    ideal-genom run --config pipeline.yaml

**Python API:**

Old::

    from ideal_genom_qc import SampleQC

    qc = SampleQC(...)
    qc.run(...)

New::

    from ideal_genom.qc.sample_qc import SampleQC

    qc = SampleQC(...)
    qc.execute_sample_qc_pipeline(...)

For detailed migration instructions, see the :doc:`getting_started` guide.

Contributing
------------

We welcome contributions! Please see :doc:`contributing` for guidelines.

For bug reports and feature requests, please use our `GitHub Issues <https://github.com/LuisGiraldo86/IDEAL-GENOM/issues>`_.
