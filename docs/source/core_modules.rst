Core Modules
============

Core framework modules for pipeline execution, configuration, and CLI.

Pipeline Executor
-----------------

.. automodule:: ideal_genom.core.pipeline
   :members:
   :undoc-members:
   :show-inheritance:

Configuration
-------------

.. automodule:: ideal_genom.core.config
   :members:
   :undoc-members:
   :show-inheritance:

Command Line Interface
----------------------

.. automodule:: ideal_genom.core.cli
   :members:
   :undoc-members:
   :show-inheritance:

Executor
--------

.. automodule:: ideal_genom.core.executor
   :members:
   :undoc-members:
   :show-inheritance:

Environment Verification
------------------------

Checks that PLINK 1.9, PLINK 2.0, GCTA and BCFtools are available on ``PATH``
and reports their versions. ``quick_verify()`` is the fastest way to confirm a
fresh installation is usable before running a pipeline.

.. automodule:: ideal_genom.core.environment
   :members:
   :undoc-members:
   :show-inheritance:

Example Datasets
----------------

Downloads the example GWAS summary-statistics datasets used by the
``viz_notebooks/`` tutorials into ``ideal_genom/data/sumstats/``. Each function
returns the :class:`~pathlib.Path` of the downloaded file and skips the
download if the file is already present.

.. automodule:: ideal_genom.core.get_examples
   :members:
   :undoc-members:
   :show-inheritance:
