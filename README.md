# IDEAL-GENOM

[![Documentation Status](https://readthedocs.org/projects/verus-ideal-genom/badge/?version=latest)](https://verus-ideal-genom.readthedocs.io/en/latest/)
[![PyPI version](https://img.shields.io/pypi/v/ideal-genom.svg)](https://pypi.org/project/ideal-genom/)
[![MIT License](https://img.shields.io/pypi/l/ideal-genom.svg)](https://pypi.org/project/ideal-genom/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://github.com/LuisGiraldo86/IDEAL-GENOM/blob/main/Dockerfile)



**IDEAL-GENOM** is a comprehensive Python package for automated, reproducible analysis of human genotype data. Currently it has implemented three pipelines: genomic quality control (QC) for case/control studies; processing VCF files after imputation; and genome-wide association studies (GWAS). It wraps years of research at CGE Tübingen, leveraging PLINK 1.9/2.0, GCTA and bcftools and also providing rich reporting and visualizations.

## Key Features

- **Sample QC**: Automated sample-level filtering.
- **Ancestry QC**: Detection of outlier ancestry samples tailored for homogenous populations and based on 1KG data.
- **Variant QC**: Automated variant-level (SNPs) filtering.
- **VCF Processing**: Post-imputation VCF to PLINK filtering, harmonization and conversion to PLINK1.9 binaries.
- **GWAS**: Generalized Linear Model (GLM) and Generalized Linear Mixed Model (GLMM) and top-hits finding.
- **Visualization**: QC steps are complemented with high quality plots for reporting. Population structure visualization powered by dimensionality reduction algorithms such as Uniform Manifold Approximation and Projection (UMAP) and t-SNE. Moreover, it has a visualization functionalities to report GWAS' summary statistics.
- **Flexible configuration**: Modular pipeline steps whose configuration is based on YAML files.
- **CLI, Jupyter, and Docker**: Run as a command-line tool, in notebooks, or containerized
- **Reproducible**: All steps, parameters, and outputs are logged.

## Installation

You can install IDEAL-GENOM using pip:

```bash
pip install ideal-genom
```

Or clone the repository and install locally:

```bash
git clone https://github.com/LuisGiraldo86/IDEAL-GENOM.git
cd IDEAL-GENOM
pip install .
```

For Docker usage:

```bash
docker build -t ideal-genom .
docker run -it ideal-genom
```

## Installed Genomic Tools in Docker

The IDEAL-GENOM Docker image comes pre-installed with the following genomic analysis tools:

- **PLINK 1.9**: Version 20231211
- **PLINK 2.0**: Version 20240105 (AVX2 build)
- **GCTA**: Version 1.95.0 (Linux x86_64)
- **BCFtools**: Version 1.23

These tools are available in the container's PATH and can be used directly in your pipeline steps or custom scripts. Example usage inside the container:

```bash
plink --help
plink2 --help
gcta64 --help
bcftools --help
```

You can run these commands interactively by starting a shell in the container:

```bash
docker run -it ideal-genom /bin/bash
```

This ensures reproducible and ready-to-use genomic analysis workflows.
