Troubleshooting Guide
=====================

This guide helps you diagnose and resolve common issues when using IDEAL-GENOM. Issues are organized by category for easier navigation.

Installation Issues
-------------------

PLINK Not Found
^^^^^^^^^^^^^^^

**Error:** ``plink: command not found`` or ``plink2: command not found``

**Solution:**

1. **Check if PLINK is installed:**

.. code-block:: bash

    which plink
    which plink2

2. **Install PLINK if missing:**

.. code-block:: bash

    # Download PLINK 1.9
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
    unzip plink_linux_x86_64_20231211.zip
    sudo mv plink /usr/local/bin/
    
    # Download PLINK 2.0
    wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20231212.zip
    unzip plink2_linux_x86_64_20231212.zip
    sudo mv plink2 /usr/local/bin/

3. **Add to PATH if installed elsewhere:**

.. code-block:: bash

    export PATH=/path/to/plink:$PATH
    # Add to ~/.bashrc for persistence

Permission Denied Errors
^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``Permission denied`` when installing or running

**Solutions:**

.. code-block:: bash

    # Install to user directory
    pip install --user ideal-genom
    
    # Or use virtual environment
    python -m venv qc_env
    source qc_env/bin/activate
    pip install ideal-genom
    
    # Fix file permissions
    chmod +x /path/to/plink

Python Module Import Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``ModuleNotFoundError: No module named 'ideal_genom_qc'``

**Solutions:**

1. **Check installation:**

.. code-block:: bash

    pip list | grep ideal
    python -c "import ideal_genom_qc; print(ideal_genom_qc.__version__)"

2. **Reinstall if needed:**

.. code-block:: bash

    pip uninstall ideal-genom
    pip install ideal-genom

3. **Check Python environment:**

.. code-block:: bash

    which python
    which pip
    # Ensure they point to the same environment

Configuration Issues
--------------------

JSON Syntax Errors
^^^^^^^^^^^^^^^^^^^

**Error:** ``JSONDecodeError: Expecting ',' delimiter``

**Solution:** Validate your JSON files:

.. code-block:: bash

    # Check JSON syntax
    python -m json.tool configFiles/parameters.json
    python -m json.tool configFiles/paths.json
    python -m json.tool configFiles/steps.json

**Common JSON mistakes:**

- Missing commas between elements
- Trailing commas after last element
- Unescaped quotes in strings
- Comments (not allowed in JSON)

File Path Issues
^^^^^^^^^^^^^^^^

**Error:** ``FileNotFoundError: [Errno 2] No such file or directory``

**Solutions:**

1. **Use absolute paths:**

.. code-block:: json

    {
        "input_directory": "/full/path/to/inputData",
        "output_directory": "/full/path/to/outputData"
    }

2. **Check file permissions:**

.. code-block:: bash

    ls -la /path/to/your/files
    # Ensure read/write permissions

3. **Verify file existence:**

.. code-block:: bash

    # Check if input files exist
    ls inputData/mydata.bed
    ls inputData/mydata.bim  
    ls inputData/mydata.fam

Invalid Parameter Values
^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``ValueError: Parameter 'mind' must be between 0 and 1``

**Solution:** Check parameter ranges:

.. code-block:: json

    {
        "sample_qc": {
            "mind": 0.2,        // Must be 0-1
            "maf": 0.01,        // Must be 0-0.5
            "hwe": 5e-8,        // Must be > 0
            "sex_check": [0.2, 0.8]  // [female_max, male_min]
        }
    }

**Parameter validation checklist:**

- ``mind``, ``geno``: 0.0 to 1.0
- ``maf``: 0.0 to 0.5  
- ``hwe``: > 0 (p-value threshold)
- ``sex_check``: [female_threshold, male_threshold] where female < male

Data Format Issues
------------------

Corrupted PLINK Files
^^^^^^^^^^^^^^^^^^^^^

**Error:** ``Error: Invalid .bed file`` or ``Error: .fam file has wrong number of columns``

**Solutions:**

1. **Validate PLINK files:**

.. code-block:: bash

    # Check file integrity
    plink --bfile inputData/mydata --freq --out test_freq
    
    # Check file formats
    head inputData/mydata.fam  # Should have 6 columns
    head inputData/mydata.bim  # Should have 6 columns
    file inputData/mydata.bed  # Should be binary

2. **Regenerate binary files:**

.. code-block:: bash

    # From PLINK text format
    plink --file inputData/mydata --make-bed --out inputData/mydata_fixed
    
    # From VCF
    plink --vcf inputData/mydata.vcf --make-bed --out inputData/mydata

Chromosome Encoding Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``Error: Unrecognized chromosome code``

**Solution:** Standardize chromosome codes:

.. code-block:: bash

    # Create chromosome update file
    echo "23 X" > update_chr.txt
    echo "24 Y" >> update_chr.txt
    echo "25 XY" >> update_chr.txt
    echo "26 MT" >> update_chr.txt
    
    # Update chromosome codes
    plink --bfile inputData/mydata --update-chr update_chr.txt --make-bed --out inputData/mydata_fixed

Missing Phenotype Data
^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``Warning: No phenotype data available``

**Solution:** Add phenotype information:

.. code-block:: bash

    # Create phenotype file (FID, IID, phenotype)
    # 1=control, 2=case, -9=missing
    awk '{print $1, $2, "1"}' inputData/mydata.fam > phenotypes.txt
    
    # Update phenotypes
    plink --bfile inputData/mydata --pheno phenotypes.txt --make-bed --out inputData/mydata_pheno

Runtime Issues
--------------

Memory Errors
^^^^^^^^^^^^^

**Error:** ``MemoryError`` or ``Killed`` (out of memory)

**Solutions:**

1. **Reduce memory usage:**

.. code-block:: json

    {
        "sample_qc": {
            "ind_pair": [200, 50, 0.2],  // Larger LD windows
            "chunk_size": 5000           // Process in chunks
        },
        "ancestry_qc": {
            "pca": 5,                    // Fewer PCs
            "maf": 0.05                  // Higher MAF filter
        }
    }

2. **Process chromosomes separately:**

.. code-block:: bash

    # Split by chromosome
    for chr in {1..22}; do
        plink --bfile inputData/mydata --chr $chr --make-bed --out chr${chr}_data
    done

3. **Monitor memory usage:**

.. code-block:: bash

    # Check available memory
    free -h
    
    # Monitor during execution
    top -p $(pgrep -f ideal_genom_qc)

Disk Space Issues
^^^^^^^^^^^^^^^^^

**Error:** ``OSError: [Errno 28] No space left on device``

**Solutions:**

1. **Check disk space:**

.. code-block:: bash

    df -h .
    du -sh outputData/

2. **Clean temporary files:**

.. code-block:: bash

    # Remove temporary PLINK files
    find . -name "*.tmp" -delete
    find . -name "plink.log" -delete
    find . -name "*.nosex" -delete

3. **Use different output directory:**

.. code-block:: json

    {
        "output_directory": "/path/to/larger/disk/outputData"
    }

Long Runtime Issues
^^^^^^^^^^^^^^^^^^^

**Issue:** Pipeline takes much longer than expected

**Solutions:**

1. **Check system resources:**

.. code-block:: bash

    # CPU usage
    htop
    
    # I/O wait
    iostat -x 1
    
    # Check for bottlenecks
    iotop

2. **Optimize parameters:**

.. code-block:: json

    {
        "sample_qc": {
            "ind_pair": [100, 25, 0.3],  // Faster LD pruning
            "use_kingship": false        // Skip if not needed
        }
    }

3. **Use SSD storage:** Move data to faster storage if possible

Output and Results Issues
-------------------------

Missing Output Files
^^^^^^^^^^^^^^^^^^^^

**Issue:** Expected output files are not generated

**Solutions:**

1. **Check pipeline logs:**

.. code-block:: bash

    # Look for error messages
    grep -i error outputData/*.log
    grep -i warning outputData/*.log

2. **Verify step completion:**

.. code-block:: bash

    # Check if steps completed
    ls outputData/*/clean_files/
    ls outputData/*/fail_samples/

3. **Re-run specific steps:**

.. code-block:: json

    {
        "ancestry": false,  // Skip completed steps
        "sample": false,
        "variant": true,    // Re-run failed step
        "umap": true
    }

Empty or Invalid Results
^^^^^^^^^^^^^^^^^^^^^^^^

**Issue:** Output files exist but are empty or contain unexpected results

**Solutions:**

1. **Check input data quality:**

.. code-block:: bash

    # Basic statistics
    plink --bfile inputData/mydata --freq --missing --out data_check
    
    # Check sample sizes
    wc -l inputData/mydata.fam
    wc -l outputData/*/clean_files/*.fam

2. **Review QC thresholds:**

.. code-block:: bash

    # Check how many samples/variants were removed
    grep -i "removed" outputData/*.log

3. **Visualize intermediate results:**

.. code-block:: python

    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Load and plot QC metrics
    metrics = pd.read_csv("outputData/sample_qc_results/qc_metrics.txt", sep="\\t")
    metrics.hist(figsize=(12, 8))
    plt.show()

Plotting and Visualization Issues
---------------------------------

Missing Plots
^^^^^^^^^^^^^

**Issue:** QC plots are not generated

**Solutions:**

1. **Check plotting dependencies:**

.. code-block:: bash

    python -c "import matplotlib, seaborn, pandas; print('All plotting modules available')"

2. **Check output directories:**

.. code-block:: bash

    ls outputData/*/plots/
    ls outputData/*_plots/

3. **Generate plots manually:**

.. code-block:: python

    from ideal_genom_qc import UMAPplot
    
    plotter = UMAPplot(
        input_path="outputData/ancestry_results/clean_files",
        input_name="clean_data",
        output_path="outputData/manual_plots"
    )
    plotter.create_umap_plots()

Plot Display Issues
^^^^^^^^^^^^^^^^^^^

**Issue:** Plots are generated but not displaying correctly

**Solutions:**

1. **Check image formats:**

.. code-block:: bash

    file outputData/*/plots/*.png
    # Should show valid PNG files

2. **Convert formats if needed:**

.. code-block:: bash

    # Convert to different format
    for img in outputData/*/plots/*.png; do
        convert "$img" "${img%.png}.pdf"
    done

3. **Check plotting backend:**

.. code-block:: python

    import matplotlib
    print(matplotlib.get_backend())
    
    # Set non-interactive backend if needed
    matplotlib.use('Agg')

Network and Download Issues
---------------------------

Reference Data Download Failures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** ``ConnectionError`` or ``TimeoutError`` when downloading reference data

**Solutions:**

1. **Check internet connection:**

.. code-block:: bash

    ping google.com
    curl -I https://github.com

2. **Manual download:**

.. code-block:: python

    from ideal_genom_qc.get_references import FetcherReference
    
    fetcher = FetcherReference(built="38")
    fetcher.download_references(
        force_redownload=True,
        timeout=300  # Increase timeout
    )

3. **Use local reference files:**

.. code-block:: json

    {
        "high_ld_file": "/path/to/local/high-LD-regions.txt"
    }

Proxy or Firewall Issues
^^^^^^^^^^^^^^^^^^^^^^^^

**Error:** Download fails due to network restrictions

**Solutions:**

1. **Configure proxy:**

.. code-block:: bash

    export http_proxy=http://proxy.company.com:8080
    export https_proxy=https://proxy.company.com:8080

2. **Download manually:** Get reference files from the GitHub repository and place them locally

Performance Optimization
------------------------

Slow Performance Debugging
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Issue:** Pipeline runs slower than expected

**Debugging steps:**

1. **Profile system resources:**

.. code-block:: bash

    # CPU and memory usage
    htop
    
    # Disk I/O
    iotop -a
    
    # Network usage (if downloading references)
    nethogs

2. **Identify bottlenecks:**

.. code-block:: python

    import cProfile
    import ideal_genom_qc
    
    # Profile the QC pipeline
    cProfile.run('ideal_genom_qc.main()', 'profile_output.txt')

3. **Optimize based on bottleneck:**

- **CPU bound:** Use fewer PCs, larger LD windows
- **Memory bound:** Process in chunks, reduce dataset size
- **I/O bound:** Use SSD, reduce intermediate file writes
- **Network bound:** Download references once, use local files

Getting Help
------------

When to Seek Additional Help
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Contact the development team if you encounter:

- Reproducible bugs not covered in this guide
- Unexpected scientific results that need expert interpretation
- Feature requests for new functionality
- Performance issues on large datasets

How to Report Issues Effectively
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When reporting issues, please include:

1. **Complete error message** (copy-paste from terminal)
2. **Configuration files** (parameters.json, paths.json, steps.json)
3. **System information:**

.. code-block:: bash

    # System info
    uname -a
    python --version
    pip show ideal-genom
    plink --version
    plink2 --version

4. **Data characteristics:**

.. code-block:: bash

    # Dataset size
    wc -l inputData/*.fam
    wc -l inputData/*.bim

5. **Steps to reproduce** the issue

6. **Expected vs. actual behavior**

**Where to get help:**

- GitHub Issues: https://github.com/LuisGiraldo86/IDEAL-GENOM/issues
- Documentation: https://verus-ideal-genom.readthedocs.io/
- Email: Contact information in the repository

Debug Mode
^^^^^^^^^^

Enable debug logging for more detailed information:

.. code-block:: python

    import logging
    logging.basicConfig(level=logging.DEBUG)
    
    # Run your QC pipeline with debug output

Or use the command line with verbose output:

.. code-block:: bash

    python -m ideal_genom_qc --verbose \\
        --path_params config/parameters.json \\
        --file_folders config/paths.json \\
        --steps config/steps.json