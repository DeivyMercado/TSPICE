#Libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import spiceypy as spy
from datetime import datetime
from tspice.utils import et_to_utc_string

#To plot a single tidal signal over time
def plot_one_signal(et, tgp, loc, colors, label=r'V_{\text{tid}}(t)/g', units = 'cm', mean_value=False, savepath=None):

    #Convert et to UTC datetime objects
    utc_times = et_to_utc_string(et)
    ticks = np.linspace(0, len(utc_times)-1, 6, dtype=int) #For x-axis ticks

    #Plot the tidal signal over time
    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(utc_times, tgp, color=colors[0])

    #Settings
    ax.set_xlabel('Time [UTC]', fontsize=12)
    ax.set_ylabel(r'$'+label+'$ ['+units+']', fontsize=12)
    ax.set_title('Lon $= %.2f$°; Lat $= %.2f$°; Depth $= %.2f$ km'%(loc['lon'], loc['lat'], loc['depth']), fontsize=12)
    ax.set_xticks([utc_times[i] for i in ticks])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.yaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.grid(alpha=0.2)
    ax.margins(x=0)

    #Plot mean value line
    if mean_value:
        tgp_mean = tgp.mean()
        ax.hlines(tgp_mean, utc_times[0], utc_times[-1], color=colors[1], label=r'Mean$=%.2f$ %s'%(tgp_mean, units))
        ax.legend(fontsize=12, loc='upper right')

    fig.tight_layout()

    #Save figure
    if savepath:
        fig.savefig(savepath, format='pdf', bbox_inches='tight', dpi=300)

#To plot multiple tidal signals on the same plot
def plot_many_signal(et, tgps, loc, colors, y_label=r'V(t)/g', signal_labels=None, units='cm', savepath=None):

    #Convert ET to UTC datetime objects
    utc_times = et_to_utc_string(et)
    ticks = np.linspace(0, len(utc_times)-1, 6, dtype=int) #For x-axis ticks

    #Plot the tidal signals over time
    fig, ax = plt.subplots(figsize=(8,4))
    for i,signal in enumerate(tgps):
        ax.plot(utc_times, signal, color=colors[i], label=signal_labels[i] if signal_labels else None)

    #Settings
    ax.set_xlabel('Time [UTC]', fontsize=12)
    ax.set_ylabel(r'$'+y_label+'$ ['+units+']', fontsize=12)
    ax.set_title('Lon $= %.2f$°; Lat $= %.2f$°; Depth $= %.2f$ km'%(loc['lon'], loc['lat'], loc['depth']), fontsize=12)
    ax.set_xticks([utc_times[i] for i in ticks])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.yaxis.set_tick_params(direction='in', which='both', labelsize=12)
    ax.grid(alpha=0.2)
    ax.margins(x=0)
    if signal_labels:
        ax.legend(fontsize=12, loc='best')

    fig.tight_layout()

    #Save figure
    if savepath:
        fig.savefig(savepath, format='pdf', bbox_inches='tight', dpi=300)