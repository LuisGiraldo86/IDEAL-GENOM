Frequently Asked Questions
===========================

This page answers common questions about using IDEAL-GENOM. If you don't find your answer here, please check the :doc:`troubleshooting` guide or open an issue on GitHub.

General Questions
-----------------

**Q: What types of genomic data does IDEAL-GENOM support?**

A: IDEAL-GENOM primarily works with human genomic data in PLINK binary format (.bed, .bim, .fam files) and VCF format. It supports:

- SNP array data (e.g., Illumina, Affymetrix)
- Imputed genotype data (VCF from TOPMed, Michigan Imputation Server)
- Whole genome sequencing (WGS) data
- Whole exome sequencing (WES) data
- Targeted sequencing panels

The pipeline is optimized for autosomal chromosomes but can handle X and Y chromosomes with appropriate configuration.

**Q: Which genome builds are supported?**

A: IDEAL-GENOM supports both major human genome builds:

- GRCh37/hg19 (use ``build: "37"`` in YAML config)
- GRCh38/hg38 (use ``build: "38"`` in YAML config, recommended)

Reference files and LD regions are automatically adjusted based on the build you specify.

**Q: Can I use IDEAL-GENOM for non-human data?**

A: The current version is specifically designed for human genomic data. While the underlying algorithms could theoretically work with other species, the reference panels, LD regions, and ancestry databases are human-specific.

Installation and Setup
-----------------------

**Q: Why do I need both PLINK 1.9 and PLINK 2.0?**

A: Different QC steps require different PLINK versions:

- PLINK 1.9: Core QC operations, kinship analysis, basic statistics
- PLINK 2.0: Advanced features, faster VCF processing, improved performance for large datasets

Both tools complement each other and are required for the full pipeline functionality.

**Q: Can I install IDEAL-GENOM without admin privileges?**

A: Yes! You can install IDEAL-GENOM in user space:

.. code-block:: bash

    # Install to user directory
    pip install --user ideal-genom
    
    # Or use a virtual environment
    python -m venv ideal_qc_env
    source ideal_qc_env/bin/activate
    pip install ideal-genom

Just ensure PLINK tools are available in your PATH or specify their location.

**Q: How do I install PLINK without admin privileges?**

A: Download PLINK binaries and add them to your PATH:

.. code-block:: bash

    # Create local bin directory
    mkdir -p ~/local/bin
    
    # Download and install PLINK 1.9
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
    unzip plink_linux_x86_64_20231211.zip
    mv plink ~/local/bin/
    
    # Add to PATH
    echo 'export PATH=$HOME/local/bin:$PATH' >> ~/.bashrc
    source ~/.bashrc

Configuration and Parameters
----------------------------

**Q: How do I choose the right QC parameters for my study?**

A: Parameter selection depends on your study type:

.. list-table:: Parameter Guidelines by Study Type
   :header-rows: 1
   :widths: 25 25 25 25

   * - Parameter
     - Population Study
     - Case-Control
     - Rare Disease
   * - mind
     - 0.05-0.1
     - 0.1-0.2
     - 0.2-0.3
   * - geno
     - 0.05
     - 0.1
     - 0.2
   * - maf
     - 0.01-0.05
     - 0.01
     - 0.0-0.001
   * - hwe
     - 1e-6
     - 5e-8
     - 1e-10

Start with conservative values and relax them if you lose too many samples/variants.

**Q: What happens if I set parameters too strictly?**

A: Overly strict parameters can lead to:

- Excessive sample removal (>20% loss is concerning)
- Loss of rare variants important for your analysis
- Population bias if certain ethnic groups are disproportionately affected
- Reduced statistical power

Monitor the QC plots and logs to ensure reasonable filtering.

**Q: Can I run only specific QC steps?**

A: Yes! Use the ``steps.json`` configuration file:

.. code-block:: yaml

   pipeline:
     steps:
       - name: "sample_qc"
         enabled: true
       - name: "ancestry_qc"
         enabled: false  # Skip ancestry QC
       - name: "variant_qc"
         enabled: true
       - name: "umap_plot"
         enabled: false  # Skip UMAP plotting

Note that some steps depend on others (e.g., variant QC needs ancestry QC results).

Data and File Formats
----------------------

**Q: How do I convert my data to PLINK format?**

A: Common conversions:

.. code-block:: bash

    # VCF to PLINK
    plink --vcf mydata.vcf --make-bed --out mydata
    
    # 23andMe format to PLINK  
    plink --23file mydata.txt --make-bed --out mydata
    
    # PLINK text to binary
    plink --file mydata --make-bed --out mydata

**Q: What should I do if my .fam file doesn't have phenotype information?**

A: If your .fam file has missing phenotypes (all -9 or 0), you can:

1. **Add phenotype data manually:**

.. code-block:: bash

    # Create phenotype file (FID, IID, phenotype)
    # 1=control, 2=case in PLINK format
    echo "FAM1 SAMPLE1 1" > phenotypes.txt
    echo "FAM2 SAMPLE2 2" >> phenotypes.txt
    
    # Update .fam file
    plink --bfile mydata --pheno phenotypes.txt --make-bed --out mydata_pheno

2. **Run without case-control specific steps:**

.. code-block:: json

    {
        "umap_plot": {
            "case_control_marker": false
        }
    }

**Q: My data has non-standard chromosome coding. How do I fix this?**

A: PLINK expects standard chromosome codes (1-22, X, Y). Convert non-standard coding:

.. code-block:: bash

    # If using 23=X, 24=Y, 25=XY, 26=MT
    plink --bfile mydata --update-chr update_chr.txt --make-bed --out mydata_fixed

Where ``update_chr.txt`` contains mappings like:

.. code-block:: text

    23 X
    24 Y
    25 XY
    26 MT

Performance and Memory
----------------------

**Q: My analysis is running very slowly. How can I speed it up?**

A: Several optimization strategies:

1. **Use faster storage:** SSD instead of HDD
2. **Increase memory:** Add more RAM if possible
3. **Parallel processing:** Use multi-core systems
4. **Reduce data size:** Filter variants/samples beforehand
5. **Optimize parameters:** Larger LD pruning windows, fewer PCs

**Q: I'm getting "out of memory" errors. What should I do?**

A: Memory issues can be addressed by:

.. code-block:: yaml

   pipeline:
     steps:
       - name: "sample_qc"
         run_params:
           ind_pair: [200, 50, 0.2]  # Larger LD windows
           chunk_size: 5000          # Process in smaller chunks
       - name: "ancestry_qc"
         run_params:
           pca: 5        # Fewer PCs
           maf: 0.05     # Higher MAF threshold

Or process chromosomes separately:

.. code-block:: bash

    # Split by chromosome first
    for chr in {1..22}; do
        plink --bfile mydata --chr $chr --make-bed --out chr${chr}_data
    done

**Q: How much disk space do I need?**

A: Disk space requirements depend on your dataset size:

.. list-table:: Disk Space Requirements
   :header-rows: 1
   :widths: 25 25 25 25

   * - Dataset
     - Input Size
     - Temp Files
     - Total Needed
   * - Small (1K samples, 100K SNPs)
     - 100MB
     - 500MB
     - 1GB
   * - Medium (10K samples, 1M SNPs)
     - 1GB
     - 5GB
     - 10GB  
   * - Large (100K samples, 5M SNPs)
     - 10GB
     - 50GB
     - 100GB

Results and Interpretation
--------------------------

**Q: How do I interpret the QC plots?**

A: Key plots to examine:

1. **Heterozygosity plot:** Should show normal distribution; outliers indicate DNA quality issues
2. **Missing data plot:** Should show most samples with <10% missing data
3. **PCA plot:** Should show clear population clusters
4. **Kinship plot:** Should identify related individuals

**Q: What constitutes "normal" QC results?**

A: Typical expectations:

- **Sample removal:** 5-15% of samples failing QC
- **Variant removal:** 10-30% of variants failing QC  
- **Population outliers:** 1-5% of samples (depends on population)
- **Related individuals:** Variable (0-10% depending on study design)

**Q: Should I be concerned if many samples fail ancestry QC?**

A: High ancestry failure rates could indicate:

1. **Wrong reference population:** Check your ``reference_pop`` setting
2. **Population admixture:** Use more lenient thresholds or "ALL" reference
3. **Technical issues:** Check for batch effects or DNA quality problems
4. **Study design:** Expected in multi-ethnic studies

Quality Control Interpretation
------------------------------

**Q: A sample failed multiple QC steps. Should I remove it?**

A: Generally yes, especially if it failed:

- High missing data + ancestry outlier
- Sex discordance + high heterozygosity  
- Multiple relatedness flags
- Technical replicates with poor concordance

**Q: Can I recover samples that failed QC?**

A: Sometimes. Options include:

1. **Relaxing thresholds:** If loss is excessive
2. **Investigating causes:** Address systematic issues
3. **Manual review:** Check borderline cases individually
4. **Batch correction:** If batch effects are detected

**Q: How do I handle related individuals?**

A: Strategies for related samples:

1. **Remove one from each pair:** Keep higher call rate individual
2. **Family-based analysis:** Use appropriate statistical methods
3. **Clustering approach:** Remove minimal set to break all relationships
4. **Separate analysis:** Analyze related/unrelated separately

Technical Issues
----------------

**Q: The pipeline crashed with a PLINK error. What should I do?**

A: Common PLINK issues:

1. **Check file formats:** Ensure files are not corrupted
2. **Verify file paths:** Use absolute paths in configuration
3. **Check disk space:** Ensure sufficient space for temp files
4. **Update PLINK:** Use latest versions
5. **Check logs:** Look for specific error messages

**Q: Reference data download failed. How do I fix this?**

A: Reference data issues:

.. code-block:: python

    from ideal_genom_qc.get_references import FetcherReference
    
    # Manually download reference data
    fetcher = FetcherReference(built="38")
    fetcher.download_references(force_redownload=True)

Or provide your own reference files in the configuration.

**Q: Can I run IDEAL-GENOM on a cluster/HPC system?**

A: Yes! Example SLURM script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=ideal_qc
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=32G
    #SBATCH --time=24:00:00
    
    module load python/3.9
    module load plink/1.9
    
    python -m ideal_genom_qc \\
        --path_params config/parameters.json \\
        --file_folders config/paths.json \\
        --steps config/steps.json \\
        --built 38

Contributing and Support
------------------------

**Q: I found a bug. How do I report it?**

A: Please report bugs on our `GitHub Issues page <https://github.com/LuisGiraldo86/IDEAL-GENOM/issues>`_ with:

1. Complete error message
2. Configuration files used
3. System information (OS, Python version, PLINK versions)
4. Steps to reproduce the issue

**Q: Can I contribute to IDEAL-GENOM development?**

A: Absolutely! We welcome contributions:

1. **Bug fixes:** Submit pull requests for any bugs you fix
2. **New features:** Propose enhancements via GitHub issues first
3. **Documentation:** Help improve documentation and examples
4. **Testing:** Report issues with different data types/systems

See our :doc:`contributing` guide for details.

**Q: Is commercial use allowed?**

A: Yes, IDEAL-GENOM is open source under the MIT license, allowing commercial use. Please review the license terms in the repository for full details.