"""
TSPICE - Tidal Signal with Python and SPICE
A package for calculating the tidal potential using SPICE.
"""

#Version of the TSPICE package
__version__ = "0.0.1"

# Import key functions to make them available at the top level for users.
# This allows users to write:
# from tspice import load_kernels, calculate_tides
# instead of:
# from tspice.utils import load_kernels

from .tidal_signal import initialize, Body
from .internal_dynamics import InteriorModel
from . import kernels
from . import utils
from . import integration_tools

#List what should be imported with "from tspice import *"
__all__ = ['initialize', 'Body', 'InteriorModel', 'internal_dynamics', 'kernels', 'utils', 'integration_tools']