.PHONY: help install lint format test test.fast test.unit test.integration test.coverage coverage.unit coverage.integration bench docs

help:
	@echo "Common targets:"
	@echo "  install        - install runtime and dev dependencies"
	@echo "  lint           - run Ruff, Black, and isort checks"
	@echo "  format         - auto-format with Black and isort"
	@echo "  test.fast      - run unit tests and non-slow integration tests"
	@echo "  test           - run full unit + integration test matrix"
	@echo "  test.coverage  - generate coverage XML reports (unit + integration)"
	@echo "  docs           - build documentation site"
	@echo "  bench          - run benchmark helper (requires compressors)"

install:
	python -m pip install --upgrade pip
	pip install -r requirements.txt
	pip install .[dev]

lint:
	ruff check src tests
	black --check src tests
	isort --check-only src tests

format:
	black src tests
	isort src tests

test: test.unit test.integration

test.fast:
	python -m pytest -m "not (integration or slow)"
	python -m pytest -m "integration and not slow"

test.unit:
	python -m pytest -m "not integration"

test.integration:
	python -m pytest -m "integration"

test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=evolutionary_compression --cov-report=xml:unit.coverage.xml -m "not integration"

coverage.integration:
	python -m pytest --cov=evolutionary_compression --cov-report=xml:integration.coverage.xml -m "integration"

bench:
	python scripts/compare_compressors.py data/fixtures/small_orthogroup.fasta

docs:
	python -m pip install --upgrade pip
	pip install pipenv
	cd docs && pipenv install || true
	cd docs && pipenv run pip install sphinx sphinx_rtd_theme Pygments
	cd docs && pipenv run sphinx-build -b html . _build/html
