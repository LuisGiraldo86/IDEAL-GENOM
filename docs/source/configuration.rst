Configuration Guide
===================

This comprehensive guide explains the YAML-based configuration system in IDEAL-GENOM |release|. The configuration file controls all aspects of your genomic analysis pipeline, from data paths to QC thresholds.

Overview
--------

IDEAL-GENOM uses a **single YAML configuration file** that defines:

- Pipeline metadata (name, output directory)
- Analysis steps to execute (QC, GWAS, VCF processing)
- Parameters for each step (thresholds, options)
- Global settings (logging, resources, file handling)

**Benefits of YAML Configuration:**

- **Single Source of Truth**: All settings in one file
- **Hierarchical Structure**: Clear organization of related parameters
- **Variable Substitution**: Reference values dynamically (e.g., ``${base_output_dir}``)
- **Step Control**: Enable/disable steps without editing code
- **Self-Documenting**: Comments explain parameters inline

Configuration File Structure
-----------------------------

A configuration file has three main sections:

.. code-block:: yaml

    pipeline:
      # Pipeline metadata and steps
      name: "my_analysis"
      base_output_dir: "/path/to/output"
      steps:
        - name: "step_name"
          # Step configuration...
    
    settings:
      # Global settings
      logging: { ... }
      resources: { ... }
      files: { ... }

Getting Started with Configuration
-----------------------------------

**1. Start from a Template**

Copy a template from the repository:

.. code-block:: bash

    cp yaml_configs/qc_pipeline_config_template.yaml my_config.yaml

**2. Edit Required Fields**

At minimum, update these paths:

.. code-block:: yaml

    pipeline:
      base_output_dir: "/your/output/path"  # Where results will be saved
      steps:
        - name: "sample_qc"
          init_params:
            input_path: "/your/input/path"  # Where your PLINK files are
            input_name: "your_dataset"      # PLINK file prefix (without .bed/.bim/.fam)

**3. Validate Your Configuration**

.. code-block:: bash

    ideal-genom validate --config my_config.yaml

Pipeline Section
----------------

The ``pipeline`` section defines your analysis workflow.

Pipeline Metadata
^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    pipeline:
      name: "my_study_qc"           # Descriptive name for this analysis
      base_output_dir: "/data/output"  # Root directory for all outputs

**name** (string, required)
    A descriptive identifier for your pipeline. Used in logging and output organization.
    
**base_output_dir** (string, required)
    Absolute path where all pipeline outputs will be saved. Each step creates subdirectories here.

Pipeline Steps
^^^^^^^^^^^^^^

Steps are executed in the order listed:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          enabled: true
          module: "ideal_genom.qc.sample_qc"
          class: "SampleQC"
          init_params:
            # Parameters passed to class __init__
          execute_params:
            # Parameters passed to execute() method

**name** (string, required)
    Unique identifier for this step. Used for variable substitution and logging.

**enabled** (boolean, required)  
    Set to ``true`` to run this step, ``false`` to skip it.

**module** (string, required)
    Python module path containing the step's class.

**class** (string, required)
    Class name to instantiate for this step.

**init_params** (mapping, required)
    Parameters passed to the class constructor (``__init__``).

**execute_params** (mapping, optional)
    Parameters passed to the ``execute()`` method when running the step.

Variable Substitution
^^^^^^^^^^^^^^^^^^^^^

Reference values from elsewhere in the configuration:

.. code-block:: yaml

    pipeline:
      base_output_dir: "/data/output"
      steps:
        - name: "sample_qc"
          init_params:
            output_path: "${base_output_dir}"  # Expands to /data/output
        
        - name: "variant_qc"
          init_params:
            # Use output from previous step
            input_path: "${steps.sample_qc.clean_dir}"
            input_name: "${steps.sample_qc.output_name}"

**Available substitutions:**

- ``${base_output_dir}`` - Pipeline's base output directory
- ``${steps.STEP_NAME.ATTRIBUTE}`` - Attributes from previous steps
  - ``.clean_dir`` - Path to clean output files
  - ``.output_name`` - Output file prefix
  - ``.output_path`` - Output directory path

QC Pipeline Configuration
--------------------------

Sample QC Step
^^^^^^^^^^^^^^

Performs individual-level quality control:

.. code-block:: yaml

    - name: "sample_qc"
      enabled: true
      module: "ideal_genom.qc.sample_qc"
      class: "SampleQC"
      init_params:
        input_path: "/data/input"           # Directory containing PLINK files
        input_name: "mydata"                # PLINK file prefix
        output_path: "${base_output_dir}"   # Output directory
        output_name: "mydata_sampleQCed"    # Output file prefix
        high_ld_regions_file: "auto"        # LD regions file (or "auto" for built-in)
        build: "38"                         # Genome build: "37" or "38"
      execute_params:
        rename_snp: true                    # Rename SNPs to chr:pos format
        hh_to_missing: true                 # Convert homozygous haploid to missing
        use_kinship: true                   # Use kinship instead of IBD
        ind_pair: [50, 5, 0.2]              # LD pruning [window, step, r²]
        mind: 0.02                          # Max missing rate per individual
        sex_check: [0.2, 0.8]               # F coefficient [female_max, male_min]
        maf: 0.01                           # Minor allele frequency threshold
        het_deviation: 3                    # Heterozygosity SD threshold
        kinship: 0.354                      # Kinship coefficient threshold
        ibd_threshold: 0.185                # IBD threshold for duplicates

**init_params:**

- **input_path** (string): Directory containing input .bed/.bim/.fam files
- **input_name** (string): PLINK file prefix (e.g., "mydata" for mydata.bed)
- **output_path** (string): Where to save QC results
- **output_name** (string): Prefix for output files
- **high_ld_regions_file** (string): Path to high-LD regions file, or any nonexistent path
  such as "auto" — there's no literal sentinel check; a path that doesn't resolve to a file
  simply triggers automatic download of the built-in regions for the given ``build``
- **build** (string): Genome build version - "37" (GRCh37/hg19) or "38" (GRCh38/hg38)

**execute_params:**

- **rename_snp** (bool): Rename SNPs to chr:pos format for consistency
- **hh_to_missing** (bool): Convert heterozygous haploid calls to missing
- **use_kinship** (bool): Use KING kinship estimation (recommended over IBD)
- **ind_pair** (list[int]): LD pruning parameters [window_size_kb, step_size_kb, r²_threshold]
  
  - window_size: SNP window in variant count (default: 50)
  - step_size: Step size in variant count (default: 5)
  - r² threshold: Correlation threshold (default: 0.2)

- **mind** (float, 0-1): Maximum missing genotype rate per individual (default: 0.02 = 2%)
- **sex_check** (list[float]): F coefficient thresholds [female_max, male_min]
  
  - female_max: Maximum F for females (default: 0.2)
  - male_min: Minimum F for males (default: 0.8)
  - Samples outside these ranges fail sex check

- **maf** (float, 0-0.5): Minor allele frequency threshold for LD pruning
- **het_deviation** (float): Standard deviations from mean heterozygosity (default: 3)
- **kinship** (float): Kinship coefficient threshold for relatedness
  
  - 0.354: 1st degree relatives
  - 0.177: 2nd degree relatives  
  - 0.088: 3rd degree relatives

- **ibd_threshold** (float): IBD threshold for identifying duplicates/monozygotic twins

Ancestry QC Step
^^^^^^^^^^^^^^^^

Detects population structure and removes ancestry outliers:

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
        pca: 10                       # Number of principal components
        maf: 0.01                     # MAF threshold for PCA
        ref_threshold: 4              # SD threshold for reference outliers
        stu_threshold: 4              # SD threshold for study outliers
        reference_pop: "EUR"          # Expected population
        num_pcs: 10                   # PCs for ancestry assignment
        distance_metric: "infinity"   # Distance metric for outlier detection

**execute_params:**

- **ind_pair** (list[int]): LD pruning parameters for PCA variants
- **pca** (int): Number of principal components to compute
- **maf** (float): MAF threshold for variants included in PCA
- **ref_threshold** (float): Standard deviations for reference population outliers
- **stu_threshold** (float): Standard deviations for study population outliers
- **reference_pop** (string): Expected population ancestry

  - "EUR": European
  - "AFR": African
  - "AMR": Admixed American
  - "EAS": East Asian
  - "SAS": South Asian

- **num_pcs** (int): Number of PCs used for ancestry classification
- **distance_metric** (string): "euclidean", "manhattan", or "infinity" (Chebyshev)

Variant QC Step
^^^^^^^^^^^^^^^

Performs variant-level quality control:

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
        miss_data_rate: 0.02          # Max missing rate across samples
        diff_genotype_rate: 1.0e-5    # Differential missingness p-value
        geno: 0.02                    # Max missing rate per variant
        maf: 0.01                     # Minor allele frequency
        hwe: 1.0e-6                   # Hardy-Weinberg equilibrium p-value
        chr_y: 24                     # Y chromosome identifier

**execute_params:**

- **miss_data_rate** (float, 0-1): Maximum overall missing data rate threshold
- **diff_genotype_rate** (float): P-value threshold for differential missingness between cases/controls
- **geno** (float, 0-1): Maximum missing genotype rate per variant
- **maf** (float, 0-0.5): Minor allele frequency threshold

  - Standard GWAS: 0.01-0.05
  - Rare variant analysis: 0.001-0.01
  - Very strict: 0.001

- **hwe** (float, 0-1): Hardy-Weinberg equilibrium p-value threshold

  - Standard: 1e-6
  - Strict: 1e-10 (for genotyping array data)
  - Relaxed: 1e-4

- **chr_y** (int): Y chromosome identifier (23 for hg19, 24 for hg38)

Population Analysis Step
^^^^^^^^^^^^^^^^^^^^^^^^

Performs dimensionality reduction and population visualization:

.. code-block:: yaml

    - name: "dimensionality_reduction"
      enabled: true
      module: "ideal_genom.population.projection"
      class: "DimensionalityReductionPipeline"
      init_params:
        input_path: "${steps.variant_qc.clean_dir}"
        input_name: "${steps.variant_qc.output_name}"
        output_path: "${base_output_dir}"
        build: "38"
        high_ld_regions_file: "auto"
        generate_plot: true
      execute_params:
        # PCA parameters
        pca_params:
          pca: 10
        force_pca_recompute: false
        
        # UMAP parameters
        run_umap: true
        umap_params:
          n_neighbors: 15
          min_dist: 0.1
          n_components: 2
        
        # t-SNE parameters
        run_tsne: true
        tsne_params:
          perplexity: 30
        
        # Plotting options
        case_control_markers: true
        plot_format: "png"
        dpi: 600

**execute_params:**

- **pca_params** (mapping): PCA configuration

  - **pca** (int): Number of components to compute

- **force_pca_recompute** (bool): Recompute PCA even if results exist
- **run_umap** (bool): Enable UMAP analysis
- **umap_params** (mapping): UMAP configuration

  - **n_neighbors** (int): Number of neighbors (5-50, default: 15)
  - **min_dist** (float): Minimum distance (0.0-1.0, default: 0.1)
  - **n_components** (int): Output dimensions (typically 2 or 3)

- **run_tsne** (bool): Enable t-SNE analysis
- **tsne_params** (mapping): t-SNE configuration

  - **perplexity** (int): Perplexity value (5-50, default: 30)

- **case_control_markers** (bool): Color by case/control status
- **plot_format** (string): "png", "svg", or "pdf"
- **dpi** (int): Plot resolution (default: 600)

Settings Section
----------------

Global settings that apply to the entire pipeline:

Logging Settings
^^^^^^^^^^^^^^^^

.. code-block:: yaml

    settings:
      logging:
        level: "INFO"              # Logging verbosity

**level** (string): Log message detail level

- "DEBUG": Very detailed, for troubleshooting
- "INFO": Standard informational messages (recommended)
- "WARNING": Only warnings and errors
- "ERROR": Only errors

.. note::

    ``settings.logging.level`` is the only logging setting actually read
    (by the ``ideal-genom run`` CLI command, which can also be overridden
    with ``--log-level``). Keys like ``file_logging`` or
    ``console_logging`` are not implemented and have no effect if present
    in a config file.

Resource Settings
^^^^^^^^^^^^^^^^^

Thread and memory usage for PLINK/GCTA subprocess calls is always sized
automatically from host resources (see ``core/utils.py``'s
``get_optimal_threads()``/``get_available_memory()``). There is currently
no ``settings.resources`` key that overrides this — a ``max_memory`` or
``max_threads`` entry under ``settings`` is accepted by the YAML parser
but silently ignored by the pipeline executor.

File Management Settings
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    settings:
      files:
        keep_intermediate: true    # Preserve temporary files

**keep_intermediate** (bool): Keep temporary intermediate files

- ``true``: Keep all files (useful for debugging)
- ``false``: Clean up after each step (saves disk space)

.. note::

    ``keep_intermediate`` is the only key under ``settings.files`` that is
    actually consumed by ``PipelineExecutor``. ``compress_outputs`` and
    ``overwrite_existing`` are not implemented and have no effect if
    present in a config file.

Report Generation Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    settings:
      reports:
        generate_reports: true     # Generate visualization reports
        plot_format: "png"         # Plot file format

**generate_reports** (bool): Automatically generate QC plots and reports

**plot_format** (string): Output format for plots

- "png": Standard format, good quality
- "svg": Vector format, scalable
- "pdf": Publication-ready format

Advanced Configuration Patterns
--------------------------------

Conditional Step Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Skip steps based on your needs:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          enabled: true
        - name: "ancestry_qc"
          enabled: false  # Skip for homogeneous population
        - name: "variant_qc"
          enabled: true
          init_params:
            # Connect directly to sample QC
            input_path: "${steps.sample_qc.clean_dir}"

Using Pre-existing Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Resume pipeline from intermediate step:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          enabled: false  # Already completed
        - name: "variant_qc"
          enabled: true
          init_params:
            # Use existing sample QC output
            input_path: "/data/output/my_study/sample_qc/clean_files"
            input_name: "mydata_sampleQCed"

Multiple Output Directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Organize outputs by analysis type:

.. code-block:: yaml

    pipeline:
      base_output_dir: "/data/project"
      steps:
        - name: "sample_qc"
          init_params:
            output_path: "${base_output_dir}/qc_results"
        - name: "gwas_prep"
          init_params:
            output_path: "${base_output_dir}/gwas_analysis"

Parameter Tuning Guidelines
----------------------------

Sample QC Thresholds
^^^^^^^^^^^^^^^^^^^^

**For Standard Case-Control GWAS:**

- mind: 0.02 (2% missing)
- maf: 0.01 (1% MAF)
- het_deviation: 3 SD
- kinship: 0.354 (exclude 1st degree relatives)

**For Rare Variant Analysis:**

- mind: 0.01 (stricter)
- maf: 0.001 (include rare variants)
- het_deviation: 4 SD (more lenient)

**For Family-Based Studies:**

- kinship: 0.088 (allow up to 3rd degree relatives)
- Adjust sex_check if samples include children

Ancestry QC Thresholds
^^^^^^^^^^^^^^^^^^^^^^^

**For Homogeneous Populations:**

- ref_threshold: 6 SD (softer)
- stu_threshold: 6 SD (softer)
- Consider disabling ancestry QC entirely

Variant QC Thresholds
^^^^^^^^^^^^^^^^^^^^^^

**For Array-Based Data:**

- geno: 0.02 (2% missing)
- hwe: 1e-10 (very strict)
- maf: 0.01

**For Sequencing Data:**

- geno: 0.05 (more lenient)
- hwe: 1e-6 (standard)
- maf: 0.001 (include rare variants)

Common Configuration Examples
------------------------------

Minimal QC Pipeline
^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    pipeline:
      name: "minimal_qc"
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
            output_name: "mydata_clean"
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            mind: 0.02
            maf: 0.01

    settings:
      logging:
        level: "INFO"

Complete QC with Ancestry
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    pipeline:
      name: "full_qc"
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
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            mind: 0.02
            sex_check: [0.2, 0.8]
            maf: 0.01
            het_deviation: 3
            kinship: 0.354
        
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
            pca: 10
            ref_threshold: 4
            stu_threshold: 4
            reference_pop: "EUR"
        
        - name: "variant_qc"
          enabled: true
          module: "ideal_genom.qc.variant_qc"
          class: "VariantQC"
          init_params:
            input_path: "${steps.ancestry_qc.clean_dir}"
            input_name: "${steps.ancestry_qc.output_name}"
            output_path: "${base_output_dir}"
            output_name: "mydata_final"
          execute_params:
            geno: 0.02
            maf: 0.01
            hwe: 1.0e-6

Troubleshooting Configuration
------------------------------

**Configuration validation fails:**

1. Check YAML syntax (indentation, colons, quotes)
2. Verify all required fields are present
3. Ensure paths exist and are accessible
4. Check module and class names are correct

**Pipeline runs but produces no output:**

1. Verify ``enabled: true`` for desired steps
2. Check input file paths are correct
3. Review ``pipeline.log`` for errors
4. Ensure output directory is writable

**Memory errors:**

1. Set ``max_memory`` explicitly
2. Reduce ``max_threads`` to free memory
3. Process datasets in batches
4. Enable ``keep_intermediate: false`` to save space

**Variable substitution not working:**

1. Ensure correct syntax: ``${variable_name}``
2. Check referenced step names match exactly
3. Verify step order (can't reference future steps)

See Also
--------

- :doc:`getting_started` - Quick start guide
- :doc:`examples` - Complete workflow examples
- :doc:`troubleshooting` - Detailed problem-solving
- :doc:`faq` - Frequently asked questions

**Docker Paths:**

When using Docker, paths should be relative to the container's ``/data`` directory:

.. code-block:: json

    {
        "input_directory": "/data/inputData",
        "input_prefix": "mydata",
        "output_directory": "/data/outputData",
        "output_prefix": "clean_data",
        "high_ld_file": "/data/dependables/high-LD-regions.txt"
    }

Steps Configuration
-------------------

The ``steps.json`` file controls which pipeline steps to execute:

.. code-block:: json

    {
        "ancestry": true,
        "sample": true,
        "variant": true,
        "umap": true,
        "fst": true
    }

**Step Dependencies:**

- ``sample`` → ``ancestry`` → ``variant`` → ``dim reduction`` → ``fst``
- You can skip steps, but maintain dependencies
- Results from previous steps are required for subsequent steps


Advanced Configuration
----------------------

Custom LD Regions
^^^^^^^^^^^^^^^^^^

Provide your own high-LD regions file:

.. code-block:: text

    # high-LD-regions.txt format
    1   48000000    52000000    # Chromosome, start, end
    2   85000000    100000000
    6   25000000    35000000

Performance Tuning
-------------------

**Memory Optimization:**

- Increase ``ind_pair`` window size for large datasets
- Reduce ``pca`` components if memory is limited
- Process chromosomes separately for very large datasets

**Speed Optimization:**

- Use SSD storage for temporary files
- Increase available CPU cores
- Consider splitting large datasets

**Disk Space Management:**

- Monitor intermediate file sizes
- Clean up temporary files regularly
- Use compression for archival storage

Best Practices
--------------

1. **Version Control**: Keep configuration files under version control
2. **Documentation**: Document parameter choices and rationale
3. **Validation**: Always validate results visually
4. **Backup**: Keep copies of successful configurations
5. **Testing**: Test parameter changes on small datasets first

Troubleshooting
---------------

**Common Configuration Issues:**

- **Path not found**: Check absolute paths and permissions
- **Parameter out of range**: Verify threshold values are reasonable
- **JSON syntax errors**: Validate JSON format
- **Memory errors**: Reduce dataset size or adjust parameters

See the :doc:`troubleshooting` guide for more detailed solutions.