# MultiNEAs Documentation

This directory contains the documentation source files for MultiNEAs.

## Building Documentation

To build the documentation locally:

```bash
cd docs
make html
```

The built documentation will be in `docs/_build/html/`.

## Documentation Structure

- `source/`: Source files for documentation
- `_build/`: Built HTML documentation (generated, not in version control)

## Requirements

Install documentation dependencies:

```bash
pip install -r requirements-dev.txt
```

This includes Sphinx and related tools for building the documentation.
