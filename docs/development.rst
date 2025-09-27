Development Workflow
====================

Environment Setup
-----------------

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   pip install .[dev]

Make Targets
------------

|project| mirrors ClipKIT's Makefile conventions. The most common targets are
shown below; run ``make help`` for the full list.

- ``make test.fast`` – unit tests plus non-slow integration coverage.
- ``make test`` – full test suite (unit and integration markers).
- ``make test.coverage`` – generates ``unit.coverage.xml`` and
  ``integration.coverage.xml`` ready for Codecov upload.
- ``make lint`` / ``make format`` – run Ruff, Black, and isort checks or apply
  formatting.
- ``make docs`` – build the Sphinx documentation site into ``docs/_build/html``.
- ``make bench`` – execute ``scripts/compare_compressors.py`` against a sample
  alignment.

Testing Strategy
----------------

- Unit tests live under ``tests/unit`` and avoid the ``integration`` marker.
- Integration tests (``tests/integration``) use the ``@pytest.mark.integration``
  decorator and may create temporary files.
- Slow tests can be skipped via the ``slow`` marker, keeping fast and full runs
  aligned with ClipKIT's approach.

Continuous Integration
----------------------

GitHub Actions run three jobs patterned after the ClipKIT pipeline:

- ``test-fast`` matrix across macOS and multiple Python versions.
- ``test-full`` serial job that produces coverage XML files and uploads to
  Codecov (requires ``CODECOV_TOKEN`` secret if the project is private).
- ``docs`` job that builds and deploys the Sphinx site to ``gh-pages`` using the
  default ``GITHUB_TOKEN``.

Static Analysis
---------------

- Ruff enforces style and catches common mistakes (``ruff check src tests``).
- Black and isort ensure consistent formatting.
- Optional type checks can be run with ``mypy src`` once type hints reach
  critical mass.
