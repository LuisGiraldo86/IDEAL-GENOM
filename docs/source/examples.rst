Examples
========

This page provides practical examples of using IDEAL-GENOM for different types of genomic studies. Each example includes complete YAML configuration files and step-by-step instructions.

Example 1: Basic QC Pipeline
-----------------------------

This example demonstrates a standard quality control pipeline for a case-control GWAS study.

**Study Setup:**

- 2,000 samples (1,000 cases, 1,000 controls)
- 500,000 SNPs genotyped on Illumina array
- European population
- Standard QC thresholds

**Complete Configuration (qc_basic.yaml):**

.. code-block:: yaml

    pipeline:
      name: "basic_qc_pipeline"
      base_output_dir: "/data/gwas_study/qc_output"
      
      steps:
        # Step 1: Sample QC
        - name: "sample_qc"
          enabled: true
          module: "ideal_genom.qc.sample_qc"
          class: "SampleQC"
          init_params:
            input_path: "/data/gwas_study/raw_data"
            input_name: "gwas_data"
            output_path: "${base_output_dir}/sample_qc"
            output_name: "sample_clean"
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            rename_snp: true
            hh_to_missing: true
            use_kinship: true
            ind_pair: [50, 5, 0.2]
            mind: 0.1
            sex_check: [0.2, 0.8]
            maf: 0.01
            het_deviation: 3
            kinship: 0.354
            ibd_threshold: 0.185
        
        # Step 2: Ancestry QC
        - name: "ancestry_qc"
          enabled: true
          module: "ideal_genom.qc.ancestry_qc"
          class: "AncestryQC"
          init_params:
            input_path: "${steps.sample_qc.output_path}"
            input_name: "${steps.sample_qc.output_name}"
            output_path: "${base_output_dir}/ancestry_qc"
            output_name: "ancestry_clean"
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            ind_pair: [50, 5, 0.2]
            pca: 10
            maf: 0.05
            ref_threshold: 3
            stu_threshold: 3
            reference_pop: "EUR"
            num_pcs: 10
        
        # Step 3: Variant QC
        - name: "variant_qc"
          enabled: true
          module: "ideal_genom.qc.variant_qc"
          class: "VariantQC"
          init_params:
            input_path: "${steps.ancestry_qc.output_path}"
            input_name: "${steps.ancestry_qc.output_name}"
            output_path: "${base_output_dir}/variant_qc"
            output_name: "final_clean"
            high_ld_file: "data/ld_regions_files/high-LD-regions_GRCH38.txt"
          execute_params:
            chr_y: 24
            miss_data_rate: 0.1
            diff_genotype_rate: 0.0001
            geno: 0.05
            maf: 0.01
            hwe: 0.000001
    
    settings:
      logging:
        level: "INFO"
        file_logging: true
      resources:
        max_memory: null
        max_threads: null
      files:
        keep_intermediate: true

**Execution:**

.. code-block:: bash

    # Validate configuration
    ideal-genom validate --config qc_basic.yaml
    
    # Preview pipeline steps
    ideal-genom run --config qc_basic.yaml --dry-run
    
    # Execute pipeline
    ideal-genom run --config qc_basic.yaml

**Output Structure:**

.. code-block:: text

    qc_output/
    ├── sample_qc/
    │   ├── sample_clean.bed/bim/fam
    │   ├── excluded_samples.txt
    │   └── qc_report.html
    ├── ancestry_qc/
    │   ├── ancestry_clean.bed/bim/fam
    │   ├── pca_results.txt
    │   └── ancestry_plot.png
    └── variant_qc/
        ├── final_clean.bed/bim/fam
        ├── excluded_variants.txt
        └── qc_summary.txt

Example 2: Complete GWAS Workflow
----------------------------------

This example shows a full workflow from QC through GWAS analysis using linear mixed models.

**Study Setup:**

- Post-QC dataset: 1,800 samples, 450,000 SNPs
- Qualitative trait (e.g., Parkinson's disease status)
- Account for population structure with PCA
- Control for relatedness with GRM

**Configuration (gwas_complete.yaml):**

.. code-block:: yaml

    pipeline:
      name: "complete_gwas"
      base_output_dir: "/data/gwas_study/gwas_results"
      
      steps:
        # Step 1: Preparatory analysis
        - name: "gwas_prep"
          enabled: true
          module: "ideal_genom.gwas.preparatory"
          class: "Preparatory"
          init_params:
            input_path: "/data/gwas_study/qc_output/variant_qc"
            input_name: "final_clean"
            output_path: "${base_output_dir}/prep"
            output_name: "gwas_ready"
            high_ld_file: "data/ld_regions_files/high-LD-regions_GRCH38.txt"
          execute_params:
            ind_pair: [50, 5, 0.2]
            pca: 10
            maf: 0.05
        
        # Step 2: Linear Mixed Model
        - name: "gwas_glmm"
          enabled: true
          module: "ideal_genom.gwas.gen_linear_mix_model"
          class: "GWAS_GLMM"
          init_params:
            input_path: "${steps.gwas_prep.output_path}"
            input_name: "${steps.gwas_prep.output_name}"
            output_path: "${base_output_dir}/glmm"
            output_name: "glmm_results"
          execute_params:
            maf: 0.01 
    
    settings:
      logging:
        level: "INFO"
        file_logging: true
      resources:
        max_threads: 8
        max_memory: 32000

**Execution:**

.. code-block:: bash

    # Run complete GWAS pipeline
    ideal-genom run --config gwas_complete.yaml


Example 3: VCF Post-Imputation Processing
------------------------------------------

This example demonstrates processing imputed VCF files from TOPMed or Michigan Imputation Server.

**Study Setup:**

- Imputed VCF files for chromosomes 1-22
- R² quality scores from imputation
- Convert to PLINK for downstream analysis
- GRCh38 genome build

**Configuration (vcf_process.yaml):**

.. code-block:: yaml

    pipeline:
      name: "imputed_data_processing"
      base_output_dir: "/data/imputation_study/processed"
      
      steps:
        # Step 1: Process VCF files
        - name: "process_vcf"
          enabled: true
          module: "ideal_genom.post_imputation.vcf_process"
          class: "ProcessVCF"
          init_params:
            input_path: "/data/imputation_study/imputed_vcfs"
            output_path: "${base_output_dir}/vcf"
            input_name: "placeholder"
            output_name: "imputed_filtered.vcf.gz"
          execute_params:
            password: null
            r2_threshold: 0.3
            build: "38"
            ref_genome: null
            ref_annotation: "/data/references/dbSNP156_GRCh38.vcf.gz"
            max_threads: null
        
        # Step 2: Convert to PLINK
        - name: "plink_conversion"
          enabled: true
          module: "ideal_genom.post_imputation.vcf_to_plink"
          class: "GetPLINK"
          init_params:
            input_path: "${steps.process_vcf.output_path}"
            input_name: "imputed_filtered"
            output_path: "${base_output_dir}/plink"
            output_name: "imputed_plink"
          execute_params:
            double_id: true
            for_fam_update_file: null
            threads: null
            memory: null
    
    settings:
      logging:
        level: "INFO"
        file_logging: true
      files:
        keep_intermediate: true

**Execution:**

.. code-block:: bash

    # Process imputed data
    ideal-genom run --config vcf_process.yaml

Example 5: Population Structure Analysis
-----------------------------------------

This example focuses on detailed population structure analysis with Fst statistics and projection.

**Study Setup:**

- Post-QC dataset with known population labels
- Calculate Fst statistics between populations
- Project samples onto reference PCA space

**Configuration (population_analysis.yaml):**

.. code-block:: yaml

    pipeline:
      name: "population_structure"
      base_output_dir: "/data/pop_structure/output"
      
      steps:
        # Ancestry QC with PCA
        - name: "ancestry_analysis"
          enabled: true
          module: "ideal_genom.qc.ancestry_qc"
          class: "AncestryQC"
          init_params:
            input_path: "/data/pop_structure/clean_data"
            input_name: "qc_passed"
            output_path: "${base_output_dir}/ancestry"
            output_name: "ancestry_results"
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            ind_pair: [50, 5, 0.2]
            pca: 20
            maf: 0.05
            ref_threshold: 6
            stu_threshold: 6
            reference_pop: "ALL"
            num_pcs: 20
        
        # Fst calculation
        - name: "fst_calculation"
          enabled: true
          module: "ideal_genom.population.fst_stats"
          class: "FstSummary"
          init_params:
            input_path: "${steps.ancestry_analysis.output_path}"
            input_name: "${steps.ancestry_analysis.output_name}"
            output_path: "${base_output_dir}/fst"
            population_file: "/data/pop_structure/populations.txt"
          execute_params:
            pairwise: true
            window_size: 50000
        
        # Dimensionality reduction
        - name: "dimensionality_reduction"
          enabled: true
          module: "ideal_genom.population.projection"
          class: "DimensionalityReductionPipeline"
          init_params:
            input_path: "${steps.ancestry_analysis.output_path}"
            input_name: "${steps.ancestry_analysis.output_name}"
            output_path: "${base_output_dir}/projection"
            reference_pca: "${steps.ancestry_analysis.pca_file}"
          execute_params:
            num_components: 10

**Execution:**

.. code-block:: bash

    ideal-genom run --config population_analysis.yaml

Python API Examples
-------------------

Using IDEAL-GENOM Programmatically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Example 1: Running QC Steps Individually**

.. code-block:: python

    from pathlib import Path
    from ideal_genom.qc.sample_qc import SampleQC
    from ideal_genom.qc.ancestry_qc import AncestryQC
    from ideal_genom.qc.variant_qc import VariantQC
    
    # Step 1: Sample QC
    sample_qc = SampleQC(
        input_path=Path("/data/raw_data"),
        input_name="genotype_data",
        output_path=Path("/data/output/sample_qc"),
        output_name="sample_clean",
        high_ld_regions_file=Path("/data/high_ld_regions.txt"),
        build="38"
    )

    sample_qc.execute_sample_qc_pipeline({
        "rename_snp": True,
        "hh_to_missing": True,
        "use_kinship": True,
        "ind_pair": [50, 5, 0.2],
        "mind": 0.1,
        "sex_check": [0.2, 0.8],
        "maf": 0.01,
        "het_deviation": 3,
        "kinship": 0.354,
        "ibd_threshold": 0.185
    })

    # Step 2: Ancestry QC
    ancestry_qc = AncestryQC(
        input_path=Path("/data/output/sample_qc"),
        input_name="sample_clean",
        output_path=Path("/data/output/ancestry_qc"),
        output_name="ancestry_clean",
        high_ld_regions_file=Path("/data/high_ld_regions.txt"),
        build="38"
    )
    
    ancestry_qc.execute_ancestry_qc_pipeline({
        "ind_pair": [50, 5, 0.2],
        "pca": 10,
        "maf": 0.05,
        "ref_threshold": 3,
        "stu_threshold": 3,
        "reference_pop": "EUR",
        "num_pcs": 10
    })
    
    # Step 3: Variant QC
    variant_qc = VariantQC(
        input_path=Path("/data/output/ancestry_qc"),
        input_name="ancestry_clean",
        output_path=Path("/data/output/variant_qc"),
        output_name="final_clean"
    )
    
    variant_qc.execute_variant_qc_pipeline({
        "chr_y": 24,
        "miss_data_rate": 0.1,
        "diff_genotype_rate": 0.0001,
        "geno": 0.05,
        "maf": 0.01,
        "hwe": 0.000001
    })
    
    print("QC pipeline completed successfully!")

**Example 2: Custom GWAS Analysis**

.. code-block:: python

    from pathlib import Path
    from ideal_genom.gwas.preparatory import Preparatory
    from ideal_genom.gwas.gen_linear_mix_model import GWAS_GLMM
    import pandas as pd
    
    # Prepare data for GWAS
    prep = Preparatory(
        input_path=Path("/data/qc_output/variant_qc"),
        input_name="final_clean",
        output_path=Path("/data/gwas/prep"),
        output_name="gwas_ready",
        high_ld_file=Path("data/ld_regions_files/high-LD-regions_GRCH38.txt")
    )
    
    prep.execute_preparatory_pipeline({
        "ind_pair": [50, 5, 0.2],
        "pca": 10,
        "maf": 0.05
    })
    
    # Run GLMM
    glmm = GWAS_GLMM(
        input_path=Path("/data/gwas/prep"),
        input_name="gwas_ready",
        output_path=Path("/data/gwas/results"),
        output_name="glmm_results"
    )
    
    glmm.execute_gwas_glmm_pipeline({
        "maf": 0.01,
        "pruned_file": Path("/data/gwas/prep/pruned_data")
    })
    
    # Load and inspect results
    results = pd.read_csv("/data/gwas/results/glmm_results.assoc.txt", sep='\t')
    significant = results[results['p'] < 5e-8]
    print(f"Found {len(significant)} genome-wide significant variants")

**Example 3: VCF Processing Pipeline**

.. code-block:: python

    from pathlib import Path
    from ideal_genom.post_imputation.vcf_process import ProcessVCF
    from ideal_genom.post_imputation.vcf_to_plink import GetPLINK
    
    # Process VCF files
    vcf_processor = ProcessVCF(
        input_path=Path("/data/imputed_vcfs"),
        output_path=Path("/data/processed"),
        input_name="placeholder",
        output_name="imputed_clean.vcf.gz"
    )
    
    vcf_processor.execute_process_vcf_pipeline({
        "password": None,
        "r2_threshold": 0.3,
        "build": "38",
        "ref_genome": None,
        "ref_annotation": "/data/references/dbSNP.vcf.gz",
        "max_threads": 16
    })
    
    # Convert to PLINK
    plink_converter = GetPLINK(
        input_path=Path("/data/processed"),
        input_name="imputed_clean",
        output_path=Path("/data/plink_output"),
        output_name="imputed_plink"
    )
    
    plink_converter.execute_plink_conversion_pipeline({
        "double_id": True,
        "for_fam_update_file": None,
        "threads": 8,
        "memory": 32000
    })
    
    print("VCF processing and conversion completed!")

Jupyter Notebook Examples
--------------------------

The package includes interactive Jupyter notebooks in the ``notebooks/`` directory:

**Available Notebooks:**

- ``00-1KG_phenotype.ipynb``: Builds the example dataset from 1000 Genomes data
- ``01-sample_qc.ipynb``: Interactive sample QC with live plotting
- ``02-ancestry_qc.ipynb``: Ancestry QC and outlier detection against the 1000 Genomes reference
- ``03-variant_qc.ipynb``: Variant-level quality control
- ``04-dimensionality_reduc.ipynb``: Population structure via PCA, UMAP and t-SNE
- ``05-fst_statistic.ipynb``: Fst statistics between the study cohort and reference populations

The ``viz_notebooks/`` directory contains the matching examples for the GWAS
summary-statistics visualizations, using the example datasets downloaded by
:mod:`ideal_genom.core.get_examples`:

- ``manhattan.ipynb`` / ``miami.ipynb``: Manhattan and Miami plots
- ``brisbane.ipynb``: Brisbane signal-density plot
- ``qq_plot.ipynb``: QQ plots
- ``beta_beta.ipynb``: Effect-size comparison between two studies
- ``trumpet_binary.ipynb`` / ``trumpet_quantitative.ipynb``: Trumpet plots

**Notebook Features:**

- Step-by-step explanations
- Interactive parameter tuning
- Real-time visualizations
- Result interpretation guides
- Export-ready plots

Common Patterns
---------------

Pattern 1: Sequential Pipeline Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run pipelines in sequence with proper data flow:

.. code-block:: bash

    # Step 1: QC Pipeline
    ideal-genom run --config qc_pipeline.yaml
    
    # Step 2: GWAS Pipeline (uses QC output)
    ideal-genom run --config gwas_pipeline.yaml
    
    # Step 3: Population Analysis
    ideal-genom run --config population_analysis.yaml

Pattern 2: Conditional Step Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Enable/disable steps based on needs:

.. code-block:: yaml

    steps:
      - name: "sample_qc"
        enabled: true      # Always run
      
      - name: "ancestry_qc"
        enabled: true      # Run if population structure is a concern
      
      - name: "variant_qc"
        enabled: false     # Skip if already done

Best Practices
--------------

Configuration Management
^^^^^^^^^^^^^^^^^^^^^^^^

1. **Use Templates**: Start with templates from ``yaml_configs/``
2. **Version Control**: Track your YAML configurations in git
3. **Comment Parameters**: Add comments explaining non-standard values
4. **Validate First**: Always run ``ideal-genom validate`` before execution

.. code-block:: yaml

    # Good: Well-documented configuration
    execute_params:
      maf: 0.05           # Higher MAF for small sample size
      hwe: 0.000001       # Standard threshold
      het_deviation: 4    # Lenient for diverse population

Data Organization
^^^^^^^^^^^^^^^^^

Organize your project directory:

.. code-block:: text

    project/
    ├── configs/
    │   ├── qc_pipeline.yaml
    │   ├── gwas_pipeline.yaml
    │   └── vcf_pipeline.yaml
    ├── data/
    │   ├── raw/
    │   ├── processed/
    │   └── results/
    ├── scripts/
    │   ├── run_analysis.sh
    │   └── visualize_results.py
    └── notebooks/
        └── exploratory_analysis.ipynb

Next Steps
----------

- Explore the :doc:`configuration` guide for detailed parameter explanations
- Check the :doc:`troubleshooting` guide for common issues
- Review pipeline-specific documentation:
  
  - :doc:`getting_started` - Quick start guide
  - :doc:`gwas_pipeline` - GWAS analysis
  - :doc:`vcf_pipeline` - VCF processing
