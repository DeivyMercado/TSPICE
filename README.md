# TSPICE

Tidal Signal with Python and SPICE

## Description

TSPICE, is a Python package developed my Bachelor's thesis work, "Tidal Potential: Calculations", to facilitate the calculation of the tidal potential and the planetary response. The package incorporates routines and integration schemes we described there. 

## Installation

It is already available on PyPI, so you can install it easily with:

```bash
pip install tspice
```

After this, you can start using the package after importing it:

```python
import tspice
```
Once the package is installed and imported into your script, you should run the following before you start using it:

```python
tspice.initialize()
```

This downloads and loads the necessary SPICE kernels from the NAIF website. Internally, \texttt{initialize} looks for the kernels in the package's structure, and if they don't exist, it loads them with \texttt{furnsh}. If the kernels and the structure where they are stored don't exist, this function downloads them, creates the structure in the package folder, and loads them. Thus, the second time you use TSPICE, \texttt{initialize} won't download the kernels again. This guarantees that the necessary kernels are available before any calculation.