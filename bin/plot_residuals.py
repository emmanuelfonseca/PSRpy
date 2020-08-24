#! /usr/bin/python

from PSRpy.tempo.resid2 import read_resid2
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

# define command-line arguments.
parser = argparse.ArgumentParser(description="Plot data from TEMPO output.")
parser.add_argument("resid2_files", action='store', help="Input TEMPO 'resid2.tmp' file containing TOA data.")
parser.add_argument("-i", action="store", dest="info_files", help="Auxillary 'info.tmp' containing TOA labels.")
parser.add_argument("-o", action="store", default="out", dest="save_filename", help="Name of output plot file (use with --save).")
parser.add_argument("--xlim", action="store", default=[], dest="x_limits", nargs=2, type=np.float, help="If set, use supplied values to set limits on the x-axis.")
parser.add_argument("--grid", action="store_true", dest="use_grid", help="If set, add grid lines to plot.")
parser.add_argument("--save", action="store_true", dest="save_plot", help="If set, save plot to file.")

# extract argument data.
args = parser.parse_args()
resid2_file = args.resid2_files
info_file = args.info_files
save_filename = args.save_filename
x_limits = args.x_limits
use_grid = args.use_grid
save_plot = args.save_plot

# plot stuff.
plt.style.use("classic")
toa_data, labels = read_resid2(resid2_file, info_file=info_file)

if (len(labels) != 0):
    # obtain unique labela and loop over them.
    unique_labels = list(set(labels))

    for current_label in unique_labels:
        current_idx = np.where(labels == current_label)
        current_toas = (toa_data["toas"])[current_idx]
        current_residuals = (toa_data["residuals"])[current_idx] * 1e6
        current_residuals_err = (toa_data["toa_uncertainties"])[current_idx] 

        plt.errorbar(
            current_toas,
            current_residuals,
            yerr=current_residuals_err,
            fmt='o'
        )

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

if use_grid:
    plt.grid()

if save_plot:
    pass

plt.show()
