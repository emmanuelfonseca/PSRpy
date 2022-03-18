#! /usr/bin/python

from matplotlib.font_manager import FontProperties
from PSRpy.tempo.resid2 import read_resid2
from PSRpy.parfile import Parfile
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

# define command-line arguments.
parser = argparse.ArgumentParser(description="Plot data from TEMPO output.")
parser.add_argument(
    "resid2_files", 
    action='store', 
    help="Input TEMPO 'resid2.tmp' file containing TOA data."
)

parser.add_argument(
    "-i", 
    action="store", 
    dest="info_file", 
    help="Auxillary 'info.tmp' containing TOA labels."
)

parser.add_argument(
    "-o", 
    action="store", 
    default="out", 
    dest="save_filename", 
    help="Name of output plot file (use with --save)."
)

parser.add_argument(
    "--dmx", 
    action="store", 
    default=None, 
    dest="dmx_file", 
    type=str, 
    help="name of file containing DMX data in format used by TEMPO's dmxparse.py script."
)

parser.add_argument(
    "--xlim", 
    action="store", 
    default=[], 
    dest="x_limits", 
    nargs=2, 
    type=np.float, 
    help="If set, use supplied values to set limits on the x-axis."
)

parser.add_argument(
    "--grid", 
    action="store_true", 
    dest="use_grid", 
    help="If set, add grid lines to plot."
)

parser.add_argument(
    "--legend", 
    action="store_true", 
    dest="use_legend", 
    help="If set, then generate a legend for the residuals figure/panel."
)

parser.add_argument(
    "--phase_orbit",
    action="store",
    default=None,
    dest="use_legend",
    nargs=2,
    help="If set, then generate a legend for the residuals figure/panel."
)

parser.add_argument(
    "--save", 
    action="store_true", 
    dest="save_plot", 
    help="If set, save plot to file."
)

parser.add_argument(
    "--years",
    action="store_true",
    dest="time_in_years",
    help="If set, convert times to decimal years."
)

# extract argument data.
args = parser.parse_args()
resid2_file = args.resid2_files
dmx_file = args.dmx_file
info_file = args.info_file
phase_orbit = args.phase_orbit
save_filename = args.save_filename
x_limits = args.x_limits
use_grid = args.use_grid
use_legend = args.use_legend
save_plot = args.save_plot
time_in_years = args.time_in_years

# extract TOA/info data.
toa_data, labels = read_resid2(resid2_file, info_file=info_file)
dmx_err = []
dmx_mjd = []
dmx_val = []
n_panels = 1

# if DMX file is supplied, then generate two panels.
if dmx_file is not None:

    try:
        print("reading DMX information from {}...".format(dmx_file))

        for current_line in open(dmx_file, "r"):
            if "#" not in current_line:
                elem = current_line.split()
                dmx_mjd += [np.float(elem[0])]
                dmx_err += [np.float(elem[2])]
                dmx_val += [np.float(elem[1])]
 
        # adjust value of matplotlib panels if DMX data are read.
        n_panels = 2

        # print a quick summary:
        print("... successfully read {0} DMX bins".format(len(dmx_val)))

    except:
        print("WARNING: supplied DMX file cannot be understood...")
        print("... ignoring DMX data, only creating one panel")


# plot stuff.
#plt.style.use("classic")
font = FontProperties()
font.set_name("sans-serif")
font.set_size(13)

fig, axs = plt.subplots(n_panels)
x_axis_label = "MJD"

if len(labels) != 0:
    # obtain unique labela and loop over them.
    unique_labels = list(set(labels))
    print("plotting {0} sets of TOAs...".format(len(unique_labels)))
    print("... with labels {0}".format(", ".join(unique_labels)))

    for current_label in unique_labels:
        current_idx = np.where(labels == current_label)
        current_toas = (toa_data["toas"])[current_idx]
        current_residuals = (toa_data["residuals"])[current_idx] * 1e6
        current_residuals_err = (toa_data["toa_uncertainties"])[current_idx] 

        if time_in_years:
            x_axis_label = "Year"
            current_toas = (current_toas - 53005.) / 365.25 + 2004.

        if (n_panels == 2):
            axs[0].errorbar(
                current_toas,
                current_residuals,
                yerr=current_residuals_err,
                fmt='+',
                label=current_label
            )

            axs[0].set_ylabel(r"$\mathcal{R}$ ($\mu$s)", fontproperties=font)

            # now set optional features if desired.
            if use_grid:
                axs[0].grid(linestyle="--")

            if use_legend:
                axs[0].legend(loc="upper left", ncol=3)

        else:
            axs.errorbar(
                current_toas,
                current_residuals,
                yerr=current_residuals_err,
                fmt='+',
                label=current_label
            )

            axs.set_xlabel(x_axis_label, fontproperties=font)
            axs.set_ylabel(r"$\mathcal{R}$ ($\mu$s)", fontproperties=font)

            # now set optional features if desired.
            if use_grid:
                axs.grid(linestyle="--")

            if use_legend:
                axs.legend(loc="upper left", ncol=3)

else:
    plt.errorbar(
        toa_data["toas"],
        toa_data["residuals"] * 1e6,
        yerr=toa_data["toa_uncertainties"],
        color="b",
        fmt="o"
    )

    if (len(x_limits) != 0):
        plt.xlim(x_limits[0], x_limits[1])


# if DMX data are read in, create the subplot.
if (dmx_file is not None and n_panels == 2):
    times = np.array(dmx_mjd)

    if (time_in_years):
        times = (times - 53005.) / 365.25 + 2004.

    axs[1].errorbar(
        times,
        np.array(dmx_val) * 1e3,
        yerr=np.array(dmx_err) * 1e3,
        fmt='+',
        color="k"
    )

    axs[1].get_shared_x_axes().join(axs[0], axs[1])
    axs[1].set_xlabel(x_axis_label, fontproperties=font)
    axs[1].set_ylabel(r"$\Delta$DM (10$^{-3}$ pc cm$^{-3}$)", fontproperties=font)

    if (use_grid):
        axs[1].grid(linestyle="--")

# save plot to file if desired.
if (save_plot):
    plt.savefig("{0}.pdf".format(save_filename), fmt="pdf")

plt.tight_layout()
plt.show()
