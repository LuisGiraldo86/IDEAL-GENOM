Getting Started
===============

This guide will help you get up and running with IDEAL-GENOM quickly. We'll walk through setting up your first genomic analysis pipeline step by step using the new YAML-based configuration system.

Overview
--------

IDEAL-GENOM uses a modern, flexible pipeline system:

1. **Prepare Your Data**: Ensure data is in PLINK1.9 format
2. **Generate Configuration**: Create a YAML configuration file
3. **Customize Pipeline**: Edit configuration to match your needs
4. **Validate Configuration**: Check for errors before running
5. **Execute Pipeline**: Run the analysis
6. **Review Results**: Examine outputs and visualizations

The New Configuration System
-----------------------------

IDEAL-GENOM uses a **YAML-based configuration system** (introduced in v0.2.0) that replaces the previous JSON approach. Benefits include:

- **Single File**: All settings in one place (no more separate parameters.json, paths.json, steps.json)
- **Hierarchical Structure**: Clear organization of pipeline steps and parameters
- **Variable Substitution**: Reference outputs from previous steps automatically
- **Enable/Disable Steps**: Easily control which analyses to run
- **Comments**: Built-in documentation within the config file

Quick Start: 5-Minute Tutorial
-------------------------------

**1. Get a Configuration Template**

Configuration templates are included in the repository under ``yaml_configs/``:

.. code-block:: bash

    # Clone the repository (if you haven't already)
    git clone https://github.com/LuisGiraldo86/IDEAL-GENOM.git
    cd IDEAL-GENOM
    
    # Copy the QC pipeline template
    cp yaml_configs/qc_pipeline_config_template.yaml my_qc_pipeline.yaml

Available templates:
- ``qc_pipeline_config_template.yaml`` - Complete QC pipeline
- ``gwas_config_template.yaml`` - GWAS analysis pipeline
- ``vcf_config_template.yaml`` - VCF post-imputation processing

**2. Edit the Configuration**

Open ``my_qc_pipeline.yaml`` and update the paths to match your data:

.. code-block:: yaml

    pipeline:
      name: "my_study_qc"
      base_output_dir: "/path/to/output"
      
      steps:
        - name: "sample_qc"
          enabled: true
          module: "ideal_genom.qc.sample_qc"
          class: "SampleQC"
          init_params:
            input_path: "/path/to/your/data"
            input_name: "mydata"
            output_path: "${base_output_dir}"
            output_name: "mydata_sampleQCed"
            high_ld_regions_file: "/path/to/high_ld_regions.txt"
            build: "38"

**3. Validate Your Configuration**

.. code-block:: bash

    ideal-genom validate --config my_qc_pipeline.yaml

**4. Run the Pipeline**

.. code-block:: bash

    ideal-genom run --config my_qc_pipeline.yaml

That's it! The pipeline will execute all enabled steps in order.

Step-by-Step Guide
------------------

Step 1: Prepare Your Data
Step 1: Prepare Your Data
^^^^^^^^^^^^^^^^^^^^^^^^^

IDEAL-GENOM works with PLINK1.9 binary format files:

- ``.bed``: Binary genotype data
- ``.bim``: Variant information (chromosome, position, alleles, etc.)
- ``.fam``: Sample information (family ID, individual ID, phenotype, etc.)

**Convert from VCF (if needed):**

.. code-block:: bash

    plink --vcf mydata.vcf.gz --make-bed --out mydata

**Data Requirements:**

- Genome build: GRCh37 (hg19) or GRCh38 (hg38)
- For ancestry QC: 1000 Genomes reference files (auto-downloaded if not provided)
- For high LD region filtering: high-LD-regions file (included with package)

Step 2: Create Your Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Option A: Use a Template (Recommended)**

Copy one of the provided templates from the repository:

.. code-block:: bash

    # Copy the QC pipeline template
    cp yaml_configs/qc_pipeline_config_template.yaml my_qc_pipeline.yaml
    
    # Or for GWAS analysis
    cp yaml_configs/gwas_config_template.yaml my_gwas_pipeline.yaml
    
    # Or for VCF processing
    cp yaml_configs/vcf_config_template.yaml my_vcf_pipeline.yaml

**Option B: Start from Scratch**

Create a minimal configuration file:

.. code-block:: yaml

    pipeline:
      name: "my_analysis"
      base_output_dir: "/data/output"
      
      steps:
        - name: "sample_qc"
          enabled: true
          module: "ideal_genom.qc.sample_qc"
          class: "SampleQC"
          init_params:
            input_path: "/data/input"
            input_name: "mydata"
            output_path: "${base_output_dir}"
            output_name: "mydata_sampleQCed"
            high_ld_regions_file: "auto"  # Use built-in file
            build: "38"
          execute_params:
            mind: 0.02
            sex_check: [0.2, 0.8]
            maf: 0.01
            het_deviation: 3
            kinship: 0.354

    settings:
      logging:
        level: "INFO"
      resources:
        max_memory: null  # Auto-detect
        max_threads: null  # Auto-detect

Step 3: Understanding the Configuration Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The YAML configuration has three main sections:

**Pipeline Section**

.. code-block:: yaml

    pipeline:
      name: "pipeline_name"              # Descriptive name for your analysis
      base_output_dir: "/path/to/output" # All outputs will go here
      steps:                             # List of analysis steps (in order)
        - name: "step_name"
          enabled: true                  # Set to false to skip this step
          module: "ideal_genom.module"   # Python module path
          class: "ClassName"             # Class to instantiate
          init_params:                   # Parameters passed to __init__
            # ...
          execute_params:                # Parameters passed to execute()
            # ...

**Variable Substitution**

Reference values from elsewhere in the config:

.. code-block:: yaml

    pipeline:
      base_output_dir: "/data/output"
      steps:
        - name: "sample_qc"
          init_params:
            output_path: "${base_output_dir}"  # Uses /data/output
        
        - name: "variant_qc"
          init_params:
            # Use output from previous step
            input_path: "${steps.sample_qc.clean_dir}"

**Settings Section**

.. code-block:: yaml

    settings:
      logging:
        level: "INFO"                    # DEBUG, INFO, WARNING, ERROR
        file_logging: true               # Log to file
        console_logging: true            # Log to console
      
      resources:
        max_memory: null                 # null = auto-detect (uses 2/3 available)
        max_threads: null                # null = auto-detect (uses cores - 2)
      
      files:
        keep_intermediate: true          # Keep temporary files
        compress_outputs: false          # Compress output files
        overwrite_existing: false        # Overwrite existing results

Step 4: Configure Your Pipeline Steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Sample QC** - Remove low-quality samples

.. code-block:: yaml

    - name: "sample_qc"
      enabled: true
      module: "ideal_genom.qc.sample_qc"
      class: "SampleQC"
      init_params:
        input_path: "/data/input"
        input_name: "mydata"
        output_path: "${base_output_dir}"
        output_name: "mydata_sampleQCed"
        high_ld_regions_file: "auto"
        build: "38"
      execute_params:
        rename_snp: true              # Rename SNPs to chr:pos format
        hh_to_missing: true           # Convert homozygous haploid calls to missing
        use_kinship: true             # Use kinship instead of IBD
        ind_pair: [50, 5, 0.2]        # LD pruning: window, step, r² threshold
        mind: 0.02                    # Max missing rate per individual (2%)
        sex_check: [0.2, 0.8]         # F coefficient bounds [female_max, male_min]
        maf: 0.01                     # Minor allele frequency threshold
        het_deviation: 3              # Heterozygosity SD threshold
        kinship: 0.354                # Kinship coefficient (2nd degree relatives)
        ibd_threshold: 0.185          # IBD threshold for duplicate detection

**Ancestry QC** - Detect population outliers

.. code-block:: yaml

    - name: "ancestry_qc"
      enabled: true
      module: "ideal_genom.qc.ancestry_qc"
      class: "AncestryQC"
      init_params:
        input_path: "${steps.sample_qc.clean_dir}"
        input_name: "${steps.sample_qc.output_name}"
        output_path: "${base_output_dir}"
        output_name: "mydata_ancestryQCed"
        high_ld_regions_file: "auto"
        build: "38"
      execute_params:
        ind_pair: [50, 5, 0.2]        # LD pruning for PCA
        pca: 10                       # Number of PCs to compute
        maf: 0.01                     # MAF threshold
        ref_threshold: 4              # SD threshold for reference outliers
        stu_threshold: 4              # SD threshold for study outliers
        reference_pop: "EUR"          # Expected population (EUR, AFR, AMR, EAS, SAS)
        num_pcs: 10                   # Number of PCs for ancestry assignment

**Variant QC** - Remove low-quality variants

.. code-block:: yaml

    - name: "variant_qc"
      enabled: true
      module: "ideal_genom.qc.variant_qc"
      class: "VariantQC"
      init_params:
        input_path: "${steps.ancestry_qc.clean_dir}"
        input_name: "${steps.ancestry_qc.output_name}"
        output_path: "${base_output_dir}"
        output_name: "mydata_variantQCed"
      execute_params:
        miss_data_rate: 0.02          # Max missing rate across all samples
        diff_genotype_rate: 1.0e-5    # Differential missingness p-value
        geno: 0.02                    # Max missing rate per variant
        maf: 0.01                     # Minor allele frequency threshold
        hwe: 1.0e-6                   # Hardy-Weinberg equilibrium p-value
        chr_y: 24                     # Y chromosome code (24 for hg38)

Step 5: Validate Your Configuration
Step 5: Validate Your Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before running the pipeline, validate your configuration:

.. code-block:: bash

    ideal-genom validate --config qc_pipeline.yaml

This checks for:

- File paths existence
- Required parameters
- Parameter value ranges
- Module and class availability
- Configuration syntax

**Example output:**

.. code-block:: text

    ✓ Configuration file is valid
    ✓ Pipeline 'my_study_qc' configured with 3/3 enabled steps

Step 6: Run the Pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic Execution**

.. code-block:: bash

    ideal-genom run --config qc_pipeline.yaml

**Dry Run (Preview Without Executing)**

.. code-block:: bash

    ideal-genom run --config qc_pipeline.yaml --dry-run

**Example dry-run output:**

.. code-block:: text

    ============================================================
    PIPELINE SUMMARY (DRY RUN)
    ============================================================
    Pipeline Name: my_study_qc
    Output Directory: /data/output
    Total Steps: 3
    Enabled Steps: 3
    
    Enabled Steps:
      1. sample_qc (ideal_genom.qc.sample_qc.SampleQC)
      2. ancestry_qc (ideal_genom.qc.ancestry_qc.AncestryQC)
      3. variant_qc (ideal_genom.qc.variant_qc.VariantQC)
    ============================================================

**Custom Logging Level**

.. code-block:: bash

    ideal-genom run --config qc_pipeline.yaml --log-level DEBUG

Step 7: Understanding the Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After pipeline execution, your output directory will contain:

.. code-block:: text

    /data/output/
    ├── my_study_qc/                    # Pipeline-specific directory
    │   ├── sample_qc/
    │   │   ├── clean_files/            # QC-passed data
    │   │   │   ├── mydata_sampleQCed.bed
    │   │   │   ├── mydata_sampleQCed.bim
    │   │   │   └── mydata_sampleQCed.fam
    │   │   ├── fail_samples/           # Removed samples with reasons
    │   │   │   ├── failed_mind.txt
    │   │   │   ├── failed_sexcheck.txt
    │   │   │   ├── failed_het.txt
    │   │   │   └── failed_kinship.txt
    │   │   └── plots/                  # Visualization reports
    │   │       ├── call_rate.png
    │   │       ├── heterozygosity.png
    │   │       ├── sex_check.png
    │   │       └── kinship_distribution.png
    │   ├── ancestry_qc/
    │   │   ├── clean_files/
    │   │   ├── fail_samples/
    │   │   │   └── ancestry_outliers.txt
    │   │   └── plots/
    │   │       ├── pca_all_samples.png
    │   │       ├── pca_after_qc.png
    │   │       └── scree_plot.png
    │   └── variant_qc/
    │       ├── clean_files/            # Final QC-passed variants
    │       │   ├── mydata_variantQCed.bed  # Ready for GWAS!
    │       │   ├── mydata_variantQCed.bim
    │       │   └── mydata_variantQCed.fam
    │       ├── fail_variants/
    │       │   ├── failed_geno.txt
    │       │   ├── failed_hwe.txt
    │       │   └── failed_maf.txt
    │       └── plots/
    │           ├── maf_distribution.png
    │           ├── hwe_distribution.png
    │           └── missingness.png
    └── pipeline.log                    # Complete execution log

**Key Output Files:**

- **clean_files/**: Final PLINK binary files ready for downstream analysis (GWAS, etc.)
- **fail_samples/fail_variants/**: Lists of excluded samples/variants with QC failure reasons
- **plots/**: Publication-ready visualizations for QC reporting
- **pipeline.log**: Detailed log of all operations, parameters, and results

Using the Python API
---------------------

For more control, use the Python API directly:

**Basic Example**

.. code-block:: python

    from ideal_genom.core.config import load_config
    from ideal_genom.core.pipeline import PipelineExecutor
    
    # Load configuration
    config = load_config("qc_pipeline.yaml")
    
    # Create and execute pipeline
    executor = PipelineExecutor(config)
    executor.execute()

**Advanced Example with Custom Handling**

.. code-block:: python

    from ideal_genom.core.config import load_config
    from ideal_genom.core.pipeline import PipelineExecutor
    import logging
    
    # Setup custom logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Load and modify configuration
    config = load_config("qc_pipeline.yaml")
    
    # Create executor
    executor = PipelineExecutor(config, dry_run=False)
    
    # Get pipeline summary
    summary = executor.get_pipeline_summary()
    print(f"Running pipeline: {summary['pipeline_name']}")
    print(f"Enabled steps: {summary['enabled_steps']}")
    
    # Execute
    try:
        executor.execute()
        print("✓ Pipeline completed successfully!")
    except Exception as e:
        print(f"✗ Pipeline failed: {e}")
        raise

**Using Individual Modules**

.. code-block:: python

    from ideal_genom.qc.sample_qc import SampleQC
    from pathlib import Path
    
    # Initialize Sample QC
    sample_qc = SampleQC(
        input_path=Path("/data/input"),
        input_name="mydata",
        output_path=Path("/data/output"),
        output_name="mydata_sampleQCed",
        high_ld_regions_file=Path("auto"),  # nonexistent path -> triggers built-in fetch
        build="38"
    )
    
    # Run with custom parameters
    sample_qc.execute_sample_qc_pipeline(sample_params={
        "rename_snp": True,
        "mind": 0.02,
        "sex_check": [0.2, 0.8],
        "maf": 0.01,
        "het_deviation": 3,
        "kinship": 0.354
    })
    
    # Access results
    print(f"Clean data saved to: {sample_qc.clean_dir}")

Common Workflows
----------------

**Workflow 1: Complete QC Pipeline**

.. code-block:: yaml

    pipeline:
      name: "full_qc"
      base_output_dir: "/data/output"
      steps:
        - name: "sample_qc"
          enabled: true
          # ... (sample QC config)
        - name: "ancestry_qc"
          enabled: true
          # ... (ancestry QC config)
        - name: "variant_qc"
          enabled: true
          # ... (variant QC config)

**Workflow 2: Skip Ancestry QC (Homogeneous Population)**

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          enabled: true
          # ...
        - name: "ancestry_qc"
          enabled: false  # Skip ancestry analysis
        - name: "variant_qc"
          enabled: true
          init_params:
            # Connect directly to sample QC output
            input_path: "${steps.sample_qc.clean_dir}"
            input_name: "${steps.sample_qc.output_name}"

**Workflow 3: Resume from Previous Step**

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          enabled: false  # Already completed
        - name: "ancestry_qc"
          enabled: false  # Already completed
        - name: "variant_qc"
          enabled: true
          init_params:
            # Use existing ancestry QC results
            input_path: "/data/output/my_study/ancestry_qc/clean_files"
            input_name: "mydata_ancestryQCed"

Tips and Best Practices
------------------------

**Configuration Management**

- Use descriptive pipeline names
- Comment your configuration extensively
- Keep configuration files in version control (git)
- Create separate configs for different studies/populations

**Resource Management**

- Set ``max_memory`` and ``max_threads`` to ``null`` for auto-detection
- For large datasets (>100K samples), consider increasing memory allocation
- Monitor logs for memory/performance issues

**Quality Control Thresholds**

- Standard thresholds work for most datasets
- For rare variant analysis, lower MAF thresholds (e.g., 0.001)
- For array data, stricter HWE thresholds (1e-10)
- Adjust kinship threshold based on study design (family vs. unrelated)

**File Organization**

- Use consistent naming conventions
- Keep intermediate files during initial runs (``keep_intermediate: true``)
- Enable logging to files (``file_logging: true``)
- Generate visualization reports (``generate_reports: true``)

**Debugging**

- Always validate configuration before running
- Use ``--dry-run`` to preview pipeline execution
- Set ``--log-level DEBUG`` for detailed troubleshooting
- Check fail_samples/fail_variants files to understand QC failures

**Debugging**

- Always validate configuration before running
- Use ``--dry-run`` to preview pipeline execution
- Set ``--log-level DEBUG`` for detailed troubleshooting
- Check fail_samples/fail_variants files to understand QC failures

Troubleshooting Common Issues
------------------------------

**Issue: "Module not found" error**

.. code-block:: text

    Solution: Check that the module path in your config is correct.
    Example: "ideal_genom.qc.sample_qc" not "ideal_genom_qc.sample_qc"

**Issue: "File not found" for input data**

.. code-block:: text

    Solution: Ensure paths are absolute or relative to execution directory.
    Use ${base_output_dir} for variable substitution.

**Issue: Pipeline runs but produces no output**

.. code-block:: text

    Solution: Check that steps are enabled: true in configuration.
    Verify input files exist at specified paths.

**Issue: High memory usage**

.. code-block:: text

    Solution: Set max_memory explicitly in settings.resources.
    Consider splitting large datasets or increasing available RAM.

Next Steps
----------

Now that you understand the basics:

- **Explore Examples**: See :doc:`examples` for complete workflows
- **Understand Configuration**: Read :doc:`configuration` for all parameters
- **Learn GWAS**: Check :doc:`gwas_pipeline` for association analysis
- **Process VCF Files**: See :doc:`vcf_pipeline` for post-imputation workflows
- **API Reference**: Browse module documentation for advanced usage

**Additional Resources:**

- Configuration templates: Clone the repository to access ``yaml_configs/`` directory
- Example notebooks in ``notebooks/`` directory
- :doc:`faq` for frequently asked questions
- :doc:`troubleshooting` for detailed problem-solving

**Getting Help:**

- GitHub Issues: https://github.com/LuisGiraldo86/IDEAL-GENOM/issues
- Check logs: Review ``pipeline.log`` for detailed execution information
- Community: Join discussions on the GitHub repository