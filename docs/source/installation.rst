Installation Guide
==================

System Requirements
-------------------

Before installing IDEAL-GENOM, ensure you have the following prerequisites:

**Software Dependencies**
   - Python 3.11 or higher (< 3.13)
   - PLINK 1.9 (required for QC and GWAS operations)
   - PLINK 2.0 (required for advanced QC operations)
   - GCTA (required for GLMM analysis)
   - BCFtools (required for VCF processing)

**Hardware Requirements**
   - Minimum 8GB RAM (16GB+ recommended for large datasets)
   - At least 20GB free disk space
   - Multi-core processor recommended for parallel processing

Installing External Tools
--------------------------

IDEAL-GENOM requires several genomic analysis tools. You can install them manually or use the provided Docker image which includes all dependencies.

**Installing PLINK**

.. code-block:: bash

    # Install PLINK 1.9
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
    unzip plink_linux_x86_64_20231211.zip
    sudo mv plink /usr/local/bin/

    # Install PLINK 2.0
    wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip
    unzip plink2_linux_avx2_20240105.zip
    sudo mv plink2 /usr/local/bin/

**Installing GCTA**

.. code-block:: bash

    # Linux
    wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-x86_64.zip
    unzip gcta-1.95.0-linux-x86_64.zip
    sudo mv gcta-1.95.0-linux-x86_64/gcta64 /usr/local/bin/

**Installing BCFtools**

.. code-block:: bash

    # Ubuntu/Debian
    sudo apt-get install bcftools
    
    # Or from source
    wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2
    tar -xjf bcftools-1.23.tar.bz2
    cd bcftools-1.23
    ./configure --prefix=/usr/local
    make && sudo make install

**Verify Installation:**

.. code-block:: bash

    plink --version
    plink2 --version
    gcta64 --version
    bcftools --version

For macOS and Windows installation instructions, please refer to the official documentation for each tool.

Installing IDEAL-GENOM
-----------------------

Option 1: PyPI Installation (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the stable version from PyPI:

.. code-block:: bash

    pip install ideal-genom

Option 2: Development Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the latest features and development version:

.. code-block:: bash

    git clone https://github.com/LuisGiraldo86/IDEAL-GENOM.git
    cd IDEAL-GENOM
    pip install -e .

Note: The development version may contain experimental features and should be used with caution in production environments.

Option 3: Docker Installation (Includes All Tools)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Docker image includes IDEAL-GENOM and all required genomic tools pre-installed:

.. code-block:: bash

    # Build from source
    git clone https://github.com/LuisGiraldo86/IDEAL-GENOM.git
    cd IDEAL-GENOM
    docker build -t ideal-genom .
    
    # Run the container
    docker run -it -v /path/to/your/data:/data ideal-genom bash
    
    # Inside the container, all tools are available:
    plink --version
    plink2 --version
    gcta64 --version
    bcftools --version
    ideal-genom --version

The Docker image includes:
   - PLINK 1.9 (v20231211)
   - PLINK 2.0 (v20240105, AVX2 build)
   - GCTA (v1.95.0)
   - BCFtools (v1.23)
   - IDEAL-GENOM and all Python dependencies

Verification
------------

Test your installation:

.. code-block:: bash

    # Check IDEAL-GENOM version
    ideal-genom --version
    
    # Generate a template to test functionality
    ideal-genom template --output test_config.yaml

**Python API test:**

.. code-block:: python

    import ideal_genom
    print(ideal_genom.__version__)
    
    from ideal_genom.core.config import load_config
    from ideal_genom.core.pipeline import PipelineExecutor
    print("✓ IDEAL-GENOM successfully installed")

Troubleshooting
---------------

**Common Issues:**

1. **External tools not found**: Ensure PLINK, GCTA, and BCFtools are installed and in your PATH, or use the Docker image
2. **Python version**: IDEAL-GENOM requires Python 3.11-3.12. Use ``python --version`` to check
3. **Permission errors**: Use ``pip install --user ideal-genom`` for user-only installation
4. **Import errors**: Ensure you're in the correct Python environment (virtual env or conda)

**Getting Help:**

- Check the :doc:`troubleshooting` guide for detailed solutions
- Report issues on `GitHub <https://github.com/LuisGiraldo86/IDEAL-GENOM/issues>`_
- See :doc:`faq` for frequently asked questions

Next Steps
----------

After installation, proceed to the :doc:`getting_started` guide to learn how to configure and run your first pipeline with IDEAL-GENOM.