# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath('../..'))
import sphinx_rtd_theme

project = 'IDEAL-GENOM'
copyright = '2026, Luis Giraldo González, Amabel Tenghe'
author = 'Luis Giraldo González, Amabel Tenghe'
release = '1.2.0'
version = '1.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.githubpages',
    'sphinx.ext.todo',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Autodoc configuration --------------------------------------------------
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'show-inheritance': True
}

autosummary_generate = True

# Mock imports for modules not declared as project dependencies.
# pandas/numpy/matplotlib/seaborn/scipy/umap/yaml/psutil are real
# pyproject.toml dependencies (always present after `poetry install`) and
# must NOT be mocked: sphinx's mock objects break mpl_toolkits.mplot3d's
# class-body tuple unpacking, which the visualization modules trigger via
# textalloc. pyarrow/Levenshtein aren't imported anywhere in ideal_genom/
# either directly or transitively in a way that needs mocking -- mocking
# pyarrow specifically breaks pandas' own internal optional-pyarrow
# version check (pandas expects a genuine ImportError when pyarrow is
# absent, not a mock object with a garbage __version__).
autodoc_mock_imports = [
    'sklearn',
    'bs4',
    'tqdm',
]

# -- Napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# -- Intersphinx configuration -----------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# -- Todo extension configuration --------------------------------------------
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"

html_theme_options = {
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_title = f"{project} v{version}"
html_short_title = project

# Add any paths that contain custom static files (such as style sheets)
html_css_files = []

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'fncychap': '',
    'maketitle': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'ideal-genom.tex', 'IDEAL-GENOM Documentation',
     'Luis Giraldo González, Amabel Tenghe', 'manual'),
]
