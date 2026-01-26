"""
tSPICE - Tidal Signal with Python and SPICE
A package for calculating the tidal potential using SPICE.
"""

#Version of the TSPICE package
__version__ = "0.0.1"

#Import key functions to make them available at the top level for users. This allows users to write:
#from tspice import initialize, Body, BodyResponse

#Instead of:
#from tspice.tidal_signal import initialize, Body
#from tspice.internal_dynamics import BodyResponse

from .tidal_signal import initialize, Body
from .internal_dynamics import BodyResponse
from . import kernels
from . import utils
from . import integration_tools
from .plotting import plot_one_signal, plot_many_signal

#List what should be imported with "from tspice import *"
__all__ = ['initialize', 'Body', 'BodyResponse', 'internal_dynamics', 'kernels', 'utils', 'integration_tools', 'plot_one_signal', 'plot_many_signal']