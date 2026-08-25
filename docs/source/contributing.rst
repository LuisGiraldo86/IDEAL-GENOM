Contributing Guide
==================

We welcome contributions to IDEAL-GENOM! This guide will help you get started with contributing to the project, whether you're fixing bugs, adding features, improving documentation, or helping with testing.

Getting Started
---------------

Development Setup
^^^^^^^^^^^^^^^^^

1. **Fork and clone the repository:**

.. code-block:: bash

    # Fork on GitHub, then clone your fork
    git clone https://github.com/YOUR_USERNAME/IDEAL-GENOM.git
    cd IDEAL-GENOM

2. **Set up development environment:**

.. code-block:: bash

    # Install Poetry (if not already installed)
    curl -sSL https://install.python-poetry.org | python3 -
    
    # Install dependencies
    poetry install
    
    # Activate virtual environment
    poetry shell

3. **Install development dependencies:**

.. code-block:: bash

    # Install additional development tools
    poetry install --with dev
    
    # Install pre-commit hooks
    pre-commit install

4. **Verify installation:**

.. code-block:: bash

    # Run tests
    pytest
    
    # Check code style
    black --check .
    flake8 .

Project Structure
^^^^^^^^^^^^^^^^^

Understanding the codebase structure:

.. code-block:: text

    IDEAL-GENOM/
    ├── ideal_genom_qc/          # Main package
    │   ├── __init__.py
    │   ├── SampleQC.py          # Sample quality control
    │   ├── AncestryQC.py        # Ancestry analysis
    │   ├── VariantQC.py         # Variant quality control
    │   ├── PopStructure.py      # Population structure analysis
    │   ├── UMAPplot.py          # UMAP visualization
    │   ├── Helpers.py           # Utility functions
    │   └── get_references.py    # Reference data handling
    ├── tests/                   # Test suite
    ├── docs/                    # Documentation
    ├── notebooks/               # Example notebooks
    ├── data/                    # Reference data
    └── pyproject.toml           # Project configuration

Types of Contributions
----------------------

We welcome several types of contributions:

Bug Reports
^^^^^^^^^^^

**Before submitting a bug report:**

- Check existing issues to avoid duplicates
- Test with the latest version
- Gather system information and error logs

**Bug report template:**

.. code-block:: text

    **Bug Description**
    A clear description of what the bug is.
    
    **To Reproduce**
    Steps to reproduce the behavior:
    1. Configuration used
    2. Command executed
    3. Error encountered
    
    **Expected Behavior**
    What you expected to happen.
    
    **Environment**
    - OS: [e.g., Ubuntu 20.04]
    - Python version: [e.g., 3.9.7]
    - IDEAL-GENOM version: [e.g., 0.1.0]
    - PLINK versions: [e.g., 1.9, 2.0]
    
    **Additional Context**
    - Configuration files
    - Log files
    - Sample data characteristics

Feature Requests
^^^^^^^^^^^^^^^^

**Feature request template:**

.. code-block:: text

    **Feature Description**
    A clear description of what you want to achieve.
    
    **Use Case**
    Why is this feature needed? What problem does it solve?
    
    **Proposed Solution**
    How would you like this implemented?
    
    **Alternatives Considered**
    What other solutions have you considered?
    
    **Additional Context**
    Any other context or screenshots about the feature request.

Code Contributions
^^^^^^^^^^^^^^^^^^

**Development workflow:**

1. **Create a feature branch:**

.. code-block:: bash

    git checkout -b feature/new-qc-method
    # or
    git checkout -b bugfix/fix-memory-leak

2. **Make your changes:**

- Follow the existing code style
- Add tests for new functionality
- Update documentation as needed
- Keep commits atomic and well-described

3. **Test your changes:**

.. code-block:: bash

    # Run all tests
    pytest
    
    # Test specific modules
    pytest tests/test_sample_qc.py
    
    # Run with coverage
    pytest --cov=ideal_genom_qc

4. **Check code quality:**

.. code-block:: bash

    # Format code
    black .
    
    # Check style
    flake8 .
    
    # Type checking
    mypy ideal_genom_qc/

5. **Commit and push:**

.. code-block:: bash

    git add .
    git commit -m "Add new QC method for contamination detection"
    git push origin feature/new-qc-method

6. **Create pull request:**

- Use the PR template
- Reference any related issues
- Include screenshots for UI changes
- Wait for review and address feedback

Documentation Contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Types of documentation improvements:**

- API documentation improvements
- Tutorial enhancements
- Example additions
- Typo fixes
- Translation (future)

**Documentation workflow:**

.. code-block:: bash

    # Install documentation dependencies
    poetry install --with docs
    
    # Build documentation locally
    cd docs/
    make html
    
    # Open in browser
    open build/html/index.html

Testing Contributions
^^^^^^^^^^^^^^^^^^^^^

**Help improve test coverage:**

.. code-block:: bash

    # Check current coverage
    pytest --cov=ideal_genom_qc --cov-report=html
    open htmlcov/index.html

**Types of tests needed:**

- Unit tests for individual functions
- Integration tests for complete workflows
- Performance tests for large datasets
- Cross-platform compatibility tests

Code Style Guidelines
---------------------

Python Style
^^^^^^^^^^^^

We follow PEP 8 with some modifications:

- **Line length:** 88 characters (Black default)
- **Imports:** Use `isort` for import sorting
- **Docstrings:** Use Google-style docstrings
- **Type hints:** Use type hints for public APIs

**Example function:**

.. code-block:: python

    def calculate_kinship_matrix(
        input_path: Path,
        output_path: Path,
        maf_threshold: float = 0.01,
        missing_threshold: float = 0.1
    ) -> pd.DataFrame:
        """Calculate kinship matrix for sample relatedness analysis.
        
        Args:
            input_path: Path to input PLINK files
            output_path: Path for output files
            maf_threshold: Minor allele frequency threshold
            missing_threshold: Maximum missing data rate
            
        Returns:
            DataFrame containing kinship coefficients
            
        Raises:
            FileNotFoundError: If input files don't exist
            ValueError: If thresholds are out of valid range
        """
        # Implementation here
        pass

Documentation Style
^^^^^^^^^^^^^^^^^^^

- **RestructuredText:** Use .rst format for documentation
- **Clear examples:** Include working code examples
- **Cross-references:** Link between related sections
- **Screenshots:** Include for UI elements

**Example documentation:**

.. code-block:: rst

    Sample Quality Control
    ======================
    
    The :class:`SampleQC` class performs comprehensive quality control
    on individual samples in your genomic dataset.
    
    Basic Usage
    -----------
    
    .. code-block:: python
    
        from pathlib import Path
        from ideal_genom.qc.sample_qc import SampleQC

        qc = SampleQC(
            input_path=Path("data/input"),
            input_name="mydata",
            output_path=Path("data/output"),
            output_name="clean_data",
            high_ld_regions_file=Path("data/high_ld_regions.txt")
        )

        qc.execute_sample_qc_pipeline(sample_params={...})

Git Workflow
------------

Branch Naming
^^^^^^^^^^^^^

Use descriptive branch names:

- `feature/add-contamination-detection`
- `bugfix/fix-memory-leak-in-pca`
- `docs/improve-api-documentation`
- `test/add-integration-tests`

Commit Messages
^^^^^^^^^^^^^^^

Follow conventional commit format:

.. code-block:: text

    type(scope): description
    
    [optional body]
    
    [optional footer]

**Examples:**

.. code-block:: text

    feat(ancestry): add support for custom reference populations
    
    fix(sample_qc): resolve memory leak in kinship calculation
    
    docs(api): add examples to SampleQC class documentation
    
    test(variant_qc): add unit tests for HWE calculation

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `test`: Tests
- `refactor`: Code refactoring
- `perf`: Performance improvement
- `style`: Code style changes

Pull Request Process
--------------------

PR Template
^^^^^^^^^^^

**Pull request template:**

.. code-block:: text

    ## Description
    Brief description of what this PR does.
    
    ## Type of Change
    - [ ] Bug fix
    - [ ] New feature
    - [ ] Documentation update
    - [ ] Performance improvement
    - [ ] Refactoring
    
    ## Testing
    - [ ] Tests pass locally
    - [ ] Added new tests for changes
    - [ ] Tested on sample datasets
    
    ## Documentation
    - [ ] Updated API documentation
    - [ ] Updated user documentation
    - [ ] Added/updated examples
    
    ## Checklist
    - [ ] Code follows style guidelines
    - [ ] Self-review completed
    - [ ] Commented hard-to-understand areas
    - [ ] No merge conflicts
    
    ## Related Issues
    Fixes #123
    Related to #456

Review Process
^^^^^^^^^^^^^^

**What reviewers look for:**

1. **Correctness:** Does the code do what it's supposed to do?
2. **Testing:** Are there adequate tests?
3. **Documentation:** Is the code well-documented?
4. **Style:** Does it follow project conventions?
5. **Performance:** Will it negatively impact performance?
6. **Compatibility:** Will it break existing functionality?

**Responding to feedback:**

- Address all comments
- Ask for clarification if needed
- Update tests and documentation
- Force-push updates to your branch

Release Process
---------------

Versioning
^^^^^^^^^^

We use semantic versioning (semver):

- **MAJOR:** Incompatible API changes
- **MINOR:** New functionality (backward compatible)
- **PATCH:** Bug fixes (backward compatible)

**Examples:**
- `0.1.0` → `0.1.1` (bug fix)
- `0.1.1` → `0.2.0` (new feature)
- `0.2.0` → `1.0.0` (major API change)

Changelog
^^^^^^^^^

We maintain a changelog following `Keep a Changelog <https://keepachangelog.com/>`_:

.. code-block:: text

    # Changelog
    
    ## [Unreleased]
    ### Added
    - New contamination detection method
    
    ### Fixed
    - Memory leak in PCA calculation
    
    ## [0.1.0] - 2025-01-15
    ### Added
    - Initial release
    - Sample QC functionality
    - Ancestry analysis
    - Variant QC
    - UMAP visualization

Community Guidelines
--------------------

Code of Conduct
^^^^^^^^^^^^^^^

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and constructive
- Welcome newcomers and help them learn
- Focus on what's best for the community
- Use inclusive language
- Be patient with questions and mistakes

Communication
^^^^^^^^^^^^^

**Preferred channels:**

- **GitHub Issues:** Bug reports, feature requests
- **GitHub Discussions:** General questions, ideas
- **Pull Request comments:** Code-specific discussions
- **Email:** Security issues, private matters

**Communication guidelines:**

- Be clear and concise
- Provide context and examples
- Use searchable, descriptive titles
- Follow up on conversations
- Tag relevant maintainers when needed

Recognition
-----------

Contributors will be recognized in:

- **Authors file:** Major contributors
- **Release notes:** Feature contributors
- **Documentation:** Example providers
- **GitHub:** All contributors via GitHub's contributor graph

**Types of recognition:**

- Code contributions
- Documentation improvements
- Bug reports and testing
- Community support
- Translations (future)

Getting Help
------------

**If you need help contributing:**

- Read existing issues and PRs for examples
- Start with "good first issue" labels
- Ask questions in GitHub discussions
- Join our community calls (when available)
- Reach out to maintainers directly

**Resources:**

- `GitHub Flow <https://guides.github.com/introduction/flow/>`_
- `Poetry documentation <https://python-poetry.org/docs/>`_
- `pytest documentation <https://docs.pytest.org/>`_
- `Sphinx documentation <https://www.sphinx-doc.org/>`_

Thank you for contributing to IDEAL-GENOM! 🎉