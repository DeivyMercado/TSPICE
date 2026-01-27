##################################################################
#                                                                #
# tSPICE: Tidal Signal with Python and SPICE                     #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)       #
# Authors: Jorge I. Zuluaga, Juanita A. Agudelo                  #
# Contact: jorge.zuluaga@udea.edu.co                             #
##################################################################
"""
Plotting Module.

This module provides functions for visualizing tidal signals.
"""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from tspice.utils import et_to_utc_string

def plot_one_signal(et, tgp, loc, colors, label=r'V_{\text{tid}}(t)/g', units='cm', mean_value=False, savepath=None):
    """
    Plot a single tidal signal time series.

    Parameters
    ----------
    et : np.ndarray
        Array of Ephemeris Times.
    tgp : np.ndarray
        Tidal potential signal array.
    loc : dict
        Location dictionary (for plot title).
    colors : list
        List of colors to use (e.g., ['blue', 'red']). First color for signal, second for mean line.
    label : str, optional
        Y-axis label (LaTeX formatted). Default: r'V_{\\text{tid}}(t)/g'.
    units : str, optional
        Units for the Y-axis label. Default: 'cm'.
    mean_value : bool, optional
        If True, plots a horizontal line at the mean value.
    savepath : str, optional
        File path to save the plot image. If None, the plot is not saved.
    """
    utc_times = et_to_utc_string(et)
    ticks = np.linspace(0, len(utc_times)-1, 6, dtype=int)

    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(utc_times, tgp, color=colors[0])

    ax.set_xlabel('Time [UTC]', fontsize=12)
    ax.set_ylabel(r'$'+label+'$ ['+units+']', fontsize=12)
    ax.set_title('Lon $= %.2f$°; Lat $= %.2f$°; Depth $= %.2f$ km'%(loc['lon'], loc['lat'], loc['depth']), fontsize=12)
    ax.set_xticks([utc_times[i] for i in ticks])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.yaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.grid(alpha=0.2)
    ax.margins(x=0)

    if mean_value:
        tgp_mean = tgp.mean()
        ax.hlines(tgp_mean, utc_times[0], utc_times[-1], color=colors[1], label=r'Mean$=%.2f$ %s'%(tgp_mean, units))
        ax.legend(fontsize=12, loc='upper right')

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight', dpi=300)
