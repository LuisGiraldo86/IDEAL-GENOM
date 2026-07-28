Visualization Modules
=====================

The ``ideal_genom.visualizations`` package provides functions for creating publication-ready plots for GWAS and genomic analysis.

Module Overview
---------------

.. contents:: Modules
   :local:
   :depth: 1

manhattan_type
--------------

Generate Manhattan and Miami plots for genome-wide association studies (GWAS).

**Module:** ``ideal_genom.visualizations.manhattan_type``

Features:
^^^^^^^^^

- Data processing and visualization of GWAS summary statistics
- Annotation of SNPs with gene information from various sources
- Highlighting and labeling of specific SNPs of interest
- Support for both Manhattan (single study) and Miami (two studies) plots

Key Functions:
^^^^^^^^^^^^^^

.. automodule:: ideal_genom.visualizations.manhattan_type
   :members:
   :undoc-members:
   :show-inheritance:

Usage Example:
^^^^^^^^^^^^^^

Adapted from the runnable ``viz_notebooks/manhattan.ipynb`` and ``viz_notebooks/miami.ipynb`` notebooks:

.. code-block:: python

   import pandas as pd
   from ideal_genom.visualizations.manhattan_type import manhattan_draw, miami_draw

   # Generate a Manhattan plot
   df_gwas = pd.read_csv("gwas_results.txt", sep="\t")

   manhattan_draw(
       data_df=df_gwas,
       snp_col='SNP',
       chr_col='CHR',
       pos_col='POS',
       p_col='P',
       plot_dir="./plots",
       to_highlight=pd.DataFrame(),  # no highlights
       save_name='manhattan_plot.png',
       genome_line=5e-8,
       yaxis_margin=10
   )

   # Generate a Miami plot comparing two studies
   df_top = pd.read_csv("gwas_results_top.txt", sep="\t")
   df_bottom = pd.read_csv("gwas_results_bottom.txt", sep="\t")

   miami_draw(
       df_top=df_top,
       df_bottom=df_bottom,
       snp_col='SNP',
       chr_col='CHR',
       pos_col='POS',
       p_col='P',
       plots_dir="./plots",
       legend_top='Discovery cohort',
       legend_bottom='Replication cohort',
       save_name='miami_plot.png'
   )

plots
-----

Functions for generating various plots for GWAS data analysis.

**Module:** ``ideal_genom.visualizations.plots``

Features:
^^^^^^^^^

- QQ plots for visualizing the distribution of p-values
- Beta-beta scatter plots for comparing effect sizes between studies
- Trumpet plots for visualizing power and effect sizes, for both binary and quantitative traits (single function, selected via ``mode``)

Key Functions:
^^^^^^^^^^^^^^

.. automodule:: ideal_genom.visualizations.plots
   :members:
   :undoc-members:
   :show-inheritance:

Usage Example:
^^^^^^^^^^^^^^

Adapted from the runnable ``viz_notebooks/qq_plot.ipynb``, ``viz_notebooks/beta_beta.ipynb``, ``viz_notebooks/trumpet_binary.ipynb`` and ``viz_notebooks/trumpet_quantitative.ipynb`` notebooks:

.. code-block:: python

   import pandas as pd
   from ideal_genom.visualizations.plots import qqplot_draw, beta_beta_draw, trumpet_draw

   gwas_df = pd.read_csv("gwas_results.txt", sep="\t")

   # Generate QQ plot
   qqplot_draw(
       df_gwas=gwas_df,
       plots_dir="./plots",
       pval_col='P',
       save_name='qq_plot.png',
       dpi=600
   )

   # Beta-beta plot comparing two studies
   gwas_df2 = pd.read_csv("gwas_results2.txt", sep="\t")
   beta_beta_draw(
       gwas_1=gwas_df,
       gwas_2=gwas_df2,
       snp_col='SNP',
       p_col='P',
       beta_col='BETA',
       se_col='SE',
       label_1='Discovery',
       label_2='Replication',
       plot_dir="./plots",
       save_name='beta_beta.png',
       significance=5e-8,
       annotate_coincidents=False
   )

   # Trumpet plot for a binary trait
   trumpet_draw(
       df_gwas=gwas_df,
       df_freq=None,
       plot_dir="./plots",
       snp_col='SNP',
       chr_col='CHR',
       pos_col='POS',
       maf_col='MAF',
       beta_col='BETA',
       power_ts=[0.2, 0.4, 0.6],
       n_case=1000,
       n_control=1000,
       sample_size=2000,
       sample_size_strategy='median',
       p_col='P',
       mode='binary',
       p_filter=None,
       save_name='trumpet_binary.png'
   )

   # Trumpet plot for a quantitative trait
   trumpet_draw(
       df_gwas=gwas_df,
       df_freq=None,
       plot_dir="./plots",
       snp_col='SNP',
       chr_col='CHR',
       pos_col='POS',
       maf_col='MAF',
       beta_col='BETA',
       power_ts=[0.2, 0.4, 0.6],
       sample_size=2000,
       sample_size_strategy='median',
       p_col='P',
       mode='quantitative',
       p_filter=None,
       save_name='trumpet_quantitative.png'
   )

zoom_heatmap
------------

Create zoomed heatmap visualizations of SNP associations, gene annotations, and linkage disequilibrium (LD) patterns.

**Module:** ``ideal_genom.visualizations.zoom_heatmap``

Features:
^^^^^^^^^

- Filter and annotate SNP data in a genomic region
- Calculate LD matrices using PLINK
- Generate three-panel plots with:

  1. Association plot with SNPs colored by functional consequences
  2. Gene track showing gene locations and orientations
  3. LD heatmap showing correlation patterns between SNPs

Key Functions:
^^^^^^^^^^^^^^

.. automodule:: ideal_genom.visualizations.zoom_heatmap
   :members:
   :undoc-members:
   :show-inheritance:

Usage Example:
^^^^^^^^^^^^^^

.. code-block:: python

   import pandas as pd
   from ideal_genom.visualizations.zoom_heatmap import draw_zoomed_heatmap

   # Load GWAS summary statistics
   sumstats = pd.read_csv("gwas_results.txt", sep="\t")

   # Create zoom heatmap around a lead SNP
   draw_zoomed_heatmap(
       data_df=sumstats,
       lead_snp='rs12345',
       snp_col='SNP',
       p_col='P',
       pos_col='POS',
       chr_col='CHR',
       output_folder="./plots",
       bfile_folder="data/genotypes",
       bfile_name="mydata",  # PLINK prefix, without .bed/.bim/.fam
       pval_threshold=5e-6,
       radius=500000,  # 500kb window
       build='38',
       anno_source='ensembl',
       extension='pdf'
   )

Notes
-----

**Dependencies:**
   - matplotlib
   - seaborn
   - pandas
   - numpy
   - textalloc (for label positioning)
   - pyensembl (for gene annotations)
   - PLINK 1.9 or 2.0 (for LD calculations)

**Annotation Sources:**
   Gene annotation is available from:

   - **Ensembl**: Via REST API or local GTF files
   - **RefSeq**: Via local GTF files

**Genome Builds:**
   Supported genome builds are GRCh37/hg19 ('37') and GRCh38/hg38 ('38')

**Output Formats:**
   Each function has its own default output format (e.g. ``manhattan_draw`` saves
   PNG, ``trumpet_draw`` saves PDF, ``beta_beta_draw``/``qqplot_draw`` save JPEG) —
   pass ``save_name``/``extension`` with the desired file extension to override it.
   All plots are publication-ready with customizable DPI.

See Also
--------

- :doc:`gwas_modules` - GWAS analysis modules that generate data for visualization
- :doc:`Helpers` - Annotation utilities used by visualization functions
- :doc:`api_overview` - Complete API reference
