"""
TSpy - Tidal Signal with Python
A package for calculating the tidal potential using SPICE.
"""

#Version of the TSpy package
__version__ = "0.0.1"

# Import key functions to make them available at the top level for users.
# This allows users to write:
# from tspy import load_kernels, calculate_tides
# instead of:
# from tspy.utils import load_kernels

from .tidal_signal import initialize, MainBody
from . import kernels
from . import utils

#List what should be imported with "from tspy import *"
__all__ = ['initialize', 'MainBody', 'kernels', 'utils']