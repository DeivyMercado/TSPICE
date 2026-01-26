##################################################################
# tSPICE Makefile
##################################################################

.PHONY: help install install-dev test verify clean build docs push release
RELMODE=release
PYTHON ?= python3
COMMIT_MSG ?= chore: sync tracked changes

help:
	@echo "tSPICE Development Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  install      - Install the package"
	@echo "  install-dev  - Install package in development mode with dev dependencies"
	@echo "  test         - Run tests with pytest"
	@echo "  verify       - Verify package installation"
	@echo "  clean        - Remove build artifacts and cache files"
	@echo "  build        - Build distribution packages"
	@echo ""
	@echo "  docs         - Build documentation (installs docs requirements)"
	@echo "  push         - Commit (all files) and push current branch"
	@echo "  release      - Release a new version (usage: make release RELMODE=release VERSION=x.y.z)"
	@echo "  create_env   - Create local dev environment (.tspice) and contrib directory"

create_env:
	@echo "Creating local development environment..."
	@mkdir -p contrib
	@test -d .tspice || $(PYTHON) -m venv .tspice
	@echo "Installing dependencies..."
	@. .tspice/bin/activate && pip install --upgrade pip
	@. .tspice/bin/activate && pip install -e .
	@if [ -f requirements-dev.txt ]; then . .tspice/bin/activate && pip install -r requirements-dev.txt; fi
	@echo "______________________________________________________________________"
	@echo "Environment setup complete."
	@echo "To activate the environment, run:"
	@echo "source .tspice/bin/activate"

install:
	pip install .

install-dev:
	pip install -e .
	if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

test:
	pytest

verify:
	@chmod +x bin/verify_installation.py
	@$(PYTHON) bin/verify_installation.py

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf src/*.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	find . -type d -name "tspice_data" -exec rm -rf {} +

build: clean
	$(PYTHON) -m build

docs:
	$(PYTHON) -m pip install -r docs/requirements.txt
	cd docs && $(PYTHON) -m sphinx.cmd.build -M html "." "_build"

push:
	@echo "Committing tracked changes (if any)..."
	@if ! git diff --quiet || ! git diff --cached --quiet || [ -n "$$(git status --porcelain)" ]; then \
		git add . && \
		git commit -m "$(COMMIT_MSG)"; \
	else \
		echo "Working tree is clean (tracked files); nothing to commit."; \
	fi
	@echo "Pushing current branch..."
	@git push -u origin HEAD

# Example: make release RELMODE=release VERSION=0.2.0.2
release: clean docs push
	@test -n "$(VERSION)" || (echo "ERROR: VERSION is required. Example: make release RELMODE=release VERSION=0.2.0" && exit 1)
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)