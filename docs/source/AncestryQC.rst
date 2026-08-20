Ancestry QC Module
==================

The Ancestry QC module performs population structure analysis and ancestry-based quality control.

Main Class
----------

.. autoclass:: ideal_genom.qc.ancestry_qc.AncestryQC
   :members:
   :undoc-members:
   :show-inheritance:

Reference Harmonization
-----------------------

``ReferenceGenomicMerger`` harmonizes the study data with the 1000 Genomes
reference panel before any ancestry analysis: it filters strand-ambiguous SNPs
and fixes chromosome, position and allele-flip mismatches, then merges the two
datasets.

.. autoclass:: ideal_genom.qc.ancestry_qc.ReferenceGenomicMerger
   :members:
   :undoc-members:
   :show-inheritance:

Outlier Detection
-----------------

``GenomicOutlierAnalyzer`` runs PCA on the merged dataset and flags ancestry
outliers by standardized distance from the reference and study centroids. A
sample must exceed **both** thresholds to be flagged, which keeps the method
conservative for homogeneous study populations.

.. autoclass:: ideal_genom.qc.ancestry_qc.GenomicOutlierAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Supporting Classes
------------------

.. autoclass:: ideal_genom.qc.ancestry_qc.AncestryQCReport
   :members:
   :undoc-members:
   :show-inheritance:
