API Overview
============

This section provides complete API documentation for all IDEAL-GENOM modules.

Module Organization
-------------------

IDEAL-GENOM is organized into several functional modules:

**Quality Control Modules** (``ideal_genom.qc``)
   - :doc:`SampleQC` - Sample-level quality control
   - :doc:`AncestryQC` - Population structure and ancestry analysis
   - :doc:`VariantQC` - Variant-level quality control

**GWAS Modules** (``ideal_genom.gwas``)
   - :doc:`gwas_modules` - Association analysis (GLM, GLMM, preparatory steps)

**VCF Processing Modules** (``ideal_genom.post_imputation``)
   - :doc:`vcf_modules` - VCF processing and PLINK conversion

**Population Modules** (``ideal_genom.population``)
   - :doc:`PopStructure` - FST statistics and PCA projection

**Visualization Modules** (``ideal_genom.visualizations``)
   - :doc:`visualization_modules` - Manhattan, Miami, Brisbane, QQ, beta-beta, trumpet plots, and zoom heatmaps

**Core Modules** (``ideal_genom.core``)
   - :doc:`core_modules` - Pipeline framework, configuration, CLI, environment verification, example datasets

**Utility Modules**
   - :doc:`Helpers` - Helper functions, annotations, Ensembl REST clients, power calculations
   - :doc:`get_references` - Reference data management

Quick Reference
---------------

Common Imports
^^^^^^^^^^^^^^

.. code-block:: python

   # QC modules
   from ideal_genom.qc.sample_qc import SampleQC
   from ideal_genom.qc.ancestry_qc import AncestryQC
   from ideal_genom.qc.variant_qc import VariantQC
   
   # GWAS modules
   from ideal_genom.gwas.preparatory import Preparatory
   from ideal_genom.gwas.gen_linear_model import GWAS_GLM
   from ideal_genom.gwas.gen_linear_mix_model import GWAS_GLMM
   
   # VCF processing
   from ideal_genom.post_imputation.vcf_process import ProcessVCF
   from ideal_genom.post_imputation.vcf_to_plink import GetPLINK
   
   # Population analysis
   from ideal_genom.population.fst_stats import FstSummary
   from ideal_genom.population.projection import PCAReduction, UMAPReduction, TSNEReduction
   from ideal_genom.population.projection import DimensionalityReductionPipeline
   
   # Visualization modules
   from ideal_genom.visualizations.manhattan_type import manhattan_draw, miami_draw, brisbane_draw
   from ideal_genom.visualizations.plots import qqplot_draw, beta_beta_draw, trumpet_draw
   from ideal_genom.visualizations.zoom_heatmap import draw_zoomed_heatmap
   
   # Core framework
   from ideal_genom.core.config import load_config
   from ideal_genom.core.pipeline import PipelineExecutor

Basic Usage Pattern
^^^^^^^^^^^^^^^^^^^

All analysis classes follow a consistent pattern:

.. code-block:: python

   from pathlib import Path
   from ideal_genom.qc.sample_qc import SampleQC

   # 1. Initialize the class with input/output paths
   qc = SampleQC(
       input_path=Path("data/input"),
       input_name="study_data",
       output_path=Path("data/output"),
       output_name="qc_clean",
       high_ld_regions_file=Path("data/high_ld_regions.txt"),
       build="38"
   )

   # 2. Execute the analysis with parameters
   qc.execute_sample_qc_pipeline(sample_params={
       "rename_snp": True,
       "hh_to_missing": True,
       "ind_pair": [50, 5, 0.2],
       "mind": 0.1,
       "sex_check": [0.2, 0.8],
       "maf": 0.01,
       "het_deviation": 3,
       "use_kinship": True,
       "kinship": 0.354,
       "ibd_threshold": 0.185
   })

   # 3. Access results
   print(f"Output files at: {qc.output_path}")

Pipeline Integration
^^^^^^^^^^^^^^^^^^^^

Classes can be used directly via Python API or integrated into YAML pipeline configurations:

**Python API (Direct Usage):**

.. code-block:: python

   # Instantiate and run directly
   sample_qc = SampleQC(...)
   sample_qc.execute_sample_qc_pipeline(sample_params={...})

**YAML Pipeline (Declarative Configuration):**

.. code-block:: yaml

   pipeline:
     steps:
       - name: "sample_qc"
         module: "ideal_genom.qc.sample_qc"
         class: "SampleQC"
         init_params:
           input_path: "data/input"
           input_name: "study_data"
         execute_params:
           mind: 0.1
           maf: 0.01

Module Details
--------------

See individual module documentation pages for:

- Complete class references with all methods
- Parameter descriptions and valid ranges
- Input/output file formats
- Usage examples
- Implementation notes

Index
-----

.. toctree::
   :maxdepth: 1
   
   SampleQC
   AncestryQC
   VariantQC
   gwas_modules
   vcf_modules
   PopStructure
   visualization_modules
   core_modules
   Helpers
   get_references
