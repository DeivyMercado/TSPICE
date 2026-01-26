# tSPICE

[![pypi](https://img.shields.io/badge/pypi-v0.0.2-blue)](https://pypi.org/project/tspice/) [![License](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://github.com/DeivyMercado/TSPICE/blob/master/LICENSE) [![python](https://img.shields.io/badge/python-3-grey)](https://pypi.org/project/tspice/) [![Powered by SpiceyPy](https://img.shields.io/badge/Powered%20by-SpiceyPy-blue)](https://github.com/AndrewAnnex/SpiceyPy) [![arXiv](https://img.shields.io/badge/arXiv-0000.00000-orange.svg?style=flat)](https://arxiv.org/abs/0000.00000)

<p></p>
<div align="center">
  <!-- <img src="https://raw.githubusercontent.com/DeivyMercado/TSPICE/main/docs/tspice-logo-white.webp" alt="tSPICE Logo" width="600"/> -->
  <img src="https://raw.githubusercontent.com/seap-udea/tspice/main/docs/tspice-logo-white.webp" alt="tSPICE Logo" width="600"/>
</div>
<p></p>

**tSPICE** is a Python package developed for the calculation of the tidal potential and the planetary response, using SPICE kernels.

## Description

tSPICE develops a coherent pathway from the fundamentals of tidal potential theory and elasticity to practical computations of tidal signals and elastic responses. The package uses SPICE’s kernels and modular routines to compute tidal signals and integrate planetary interior models.

This work was developed as part of the Bachelor's thesis **"Planetary Tides and Elastic Response: From Theory to Calculations"** by Deivy J. Mercado R. (2026).

Key features:
-   **Tidal Potential**: Implements Tide-Generating Potentials (TGPs) in spherical-harmonic form.
-   **Elastic Response**: Solves governing elastodynamic equations for self-gravitating, elastic, transversely isotropic, spherical bodies.
-   **Integration**: Solves coupled first-order ODEs for displacement, strain, stress, and perturbing potential fields.
-   **Validation**: Reproduces tidal signals on Earth consistent with ETERNA-x and computes Love numbers in agreement with literature (e.g., using a modified PREM Earth model).

## Installation

You can install it easily with:

```bash
pip install tspice
```

## Quick Start

After installation, you can start using the package:

```python
import tspice
```

Before your first calculation, initialize the package to load SPICE kernels:

```python
tspice.initialize()
```

This function ensures that tThe second time you use tSPICE, \texttt{initialize} won't download the kernels again.f not present) and loaded.

## Examples

*(Examples to be added)*

## Authors

- **Deivy Mercado** - david231097@gmail.com
- **Jorge I. Zuluaga** - jorge.zuluaga@udea.edu.co
- **Gloria Moncayo** - gloria.moncayo@udea.edu.co

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) - see the [LICENSE](LICENSE) file for details.