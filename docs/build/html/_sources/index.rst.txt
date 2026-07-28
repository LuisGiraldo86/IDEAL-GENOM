.. IDEAL-GENOM documentation master file

IDEAL-GENOM Documentation
==========================

.. image:: https://readthedocs.org/projects/verus-ideal-genom/badge/?version=latest
   :target: https://verus-ideal-genom.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/ideal-genom.svg
   :target: https://pypi.org/project/ideal-genom/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/l/ideal-genom.svg
   :target: https://pypi.org/project/ideal-genom/
   :alt: MIT License

.. image:: https://img.shields.io/badge/docker-ready-blue.svg
   :target: https://github.com/LuisGiraldo86/IDEAL-GENOM/blob/main/Dockerfile
   :alt: Docker



**IDEAL-GENOM** is a comprehensive Python package for automated, reproducible analysis of human genotype data. It provides end-to-end pipelines for genomic quality control (QC), post-imputation VCF processing, and genome-wide association studies (GWAS). The package wraps years of research expertise from CGE Tübingen, integrating PLINK 1.9/2.0, GCTA, and BCFtools with rich reporting and visualizations.

Version: **1.2.0**

🎯 Key Features
---------------

**Comprehensive Pipelines**
   - **Genomic QC**: Sample QC, Ancestry QC, and Variant QC for case-control studies
   - **GWAS Analysis**: Generalized Linear Models (GLM) and Mixed Models (GLMM)
   - **VCF Processing**: Post-imputation filtering, normalization, and conversion to PLINK
   - **Population Structure**: FST statistics, PCA, UMAP visualization, and ancestry projection

**Advanced Analytics**
   - **Sample Quality Control**: Missingness, sex verification, heterozygosity, relatedness (kinship/IBD)
   - **Ancestry Analysis**: Population stratification detection with 1000 Genomes reference
   - **Variant Filtering**: Hardy-Weinberg equilibrium, MAF, genotype rate, differential missingness
   - **GWAS Tools**: Association testing, top-hits extraction, gene annotation (Ensembl/RefSeq)
   - **Dimensionality Reduction**: PCA and UMAP for population structure visualization

**Modern Design**
   - **YAML Configuration**: Single configuration file with clear, hierarchical structure
   - **Flexible Pipeline System**: Enable/disable steps, customize parameters per analysis
   - **Multiple Interfaces**: Command-line tool, Python API, Jupyter notebooks
   - **Docker Support**: Pre-configured container with all genomic tools installed
   - **Automated Workflows**: Pipeline executor handles dependencies and data flow
   - **Rich Reporting**: Publication-ready plots and comprehensive QC metrics

**Modern Design**
   - **YAML Configuration**: Single configuration file with clear, hierarchical structure
   - **Flexible Pipeline System**: Enable/disable steps, customize parameters per analysis
   - **Multiple Interfaces**: Command-line tool, Python API, Jupyter notebooks
   - **Docker Support**: Pre-configured container with all genomic tools installed
   - **Automated Workflows**: Pipeline executor handles dependencies and data flow
   - **Rich Reporting**: Publication-ready plots and comprehensive QC metrics

**Developer Friendly**
   - **Reproducible**: All steps, parameters, and outputs logged
   - **Extensible**: Modular architecture for adding custom analysis steps
   - **Well Documented**: Comprehensive guides, API reference, and examples
   - **Type Hints**: Full type annotations for better IDE support

Quick Start
-----------

**Installation**

.. code-block:: bash

   pip install ideal-genom

**Basic Usage**

.. code-block:: bash

   # Generate a configuration template
   ideal-genom template --output my_pipeline.yaml
   
   # Edit the configuration file to match your data
   nano my_pipeline.yaml
   
   # Validate your configuration
   ideal-genom validate --config my_pipeline.yaml
   
   # Run the pipeline
   ideal-genom run --config my_pipeline.yaml

**Python API**

.. code-block:: python

   from ideal_genom.core.config import load_config
   from ideal_genom.core.pipeline import PipelineExecutor
   
   # Load configuration
   config = load_config("my_pipeline.yaml")
   
   # Create and execute pipeline
   executor = PipelineExecutor(config)
   executor.execute()

Available Pipelines
-------------------

**QC Pipeline** - Quality control for case-control studies
   1. Sample QC: Individual-level quality control
   2. Ancestry QC: Population structure and outlier detection
   3. Variant QC: SNP-level quality control
   4. Population Visualization: UMAP/t-SNE plots

**GWAS Pipeline** - Genome-wide association analysis
   1. Preparatory: LD pruning and PCA decomposition
   2. GLM Analysis: Fixed effects association testing
   3. GLMM Analysis: Mixed model with genetic relationship matrix
   4. Annotation: Gene mapping and functional annotation

**VCF Pipeline** - Post-imputation processing
   1. VCF Processing: Filter, normalize, annotate, concatenate
   2. PLINK Conversion: Convert to PLINK binary format
   3. Quality filtering: R² threshold, multiallelic handling

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   getting_started
   configuration
   examples
   
.. toctree::
   :maxdepth: 2
   :caption: Pipelines

   qc_pipeline
   gwas_pipeline
   vcf_pipeline

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api_overview

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   faq
   troubleshooting
   contributing
   changelog

Supported Tools
---------------

IDEAL-GENOM integrates the following genomic analysis tools:

- **PLINK 1.9**: Classic PLINK for QC and association analysis
- **PLINK 2.0**: Modern version with improved performance (AVX2 optimized)
- **GCTA**: Genetic relationship matrix and mixed model analysis
- **BCFtools**: VCF manipulation and quality filtering

These tools are automatically used by the pipeline and must be installed separately or use the provided Docker image.

Citation
--------

If you use IDEAL-GENOM in your research, please cite:

.. code-block:: bibtex

   @software{ideal_genom_2026,
     title = {IDEAL-GENOM: Comprehensive Genomic Analysis Pipeline},
     author = {Giraldo González, Luis and Tenghe, Amabel},
     year = {2026},
     version = {1.2.0},
     url = {https://github.com/cge-tubingens/ideal-genom-qc}
   }

Getting Help
------------

- **Documentation**: https://ideal-genom-qc.readthedocs.io/
- **Issues**: https://github.com/cge-tubingens/cge-comrare-pipeline/issues
- **Examples**: See the :doc:`examples` page for complete workflows

License
-------

IDEAL-GENOM is released under the MIT License. See the LICENSE file in the repository for details.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

