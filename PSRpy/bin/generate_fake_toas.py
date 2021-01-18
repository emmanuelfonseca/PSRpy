#! /usr/bin/python

from PSRpy.simulate import simulate_TOAs
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(
    description="A script that uses PSRpy API to simulate TOAs that are " + \
    "consistent with the parameters of an input timing solution"
)

parser.add_argument(
    "parfile",
    type=str,
    help="a TEMPO-compatible timing solution"
)

parser.add_argument(
    "--timespan",
    default=2.0,
    dest="timespan",
    type=np.float,
    help="the timespan of the simulated TOA data set"
)

parser.add_argument(
    "--start",
    default=58515.0,
    dest="epoch_start",
    type=np.float,
    help="the starting epoch in MJD"
)

parser.add_argument(
    "--rms",
    default=1.0,
    dest="rms_residual",
    type=np.float,
    help="the desired RMS residual of simualted data set"
)

parser.add_argument(
    "--uncertainty",
    default=1.0,
    dest="toa_uncertainty",
    type=np.float,
    help="the mean uncertainty of simulated TOAs"
)

parser.add_argument(
    "--Nchannels",
    default=1,
    dest="n_chan",
    type=np.int,
    help="the number of frequency channels to simulate"
)

parser.add_argument(
    "--Nepochs",
    default=6,
    dest="n_epochs",
    type=np.int,
    help="the number of observing epochs to simulate"
)

# now retrieve data from the command line.
args = parser.parse_args()
input_parfile = args.parfile
timespan = args.timespan
epoch_start = args.epoch_start
rms_residual = args.rms_residual
toa_uncertainty = args.toa_uncertainty
n_chan = args.n_chan
n_epochs = args.n_epochs

# derive quantities from command-line values.
epoch_finish = epoch_start + (365.25 * timespan)

# now simulate data set.
print("creating TOA set for {0} epochs...".format(n_epochs))

sim = simulate_TOAs(
    input_parfile, 
    epoch_start=epoch_start, 
    epoch_finish=epoch_finish, 
    n_channels_per_epoch=n_chan, 
    n_epochs=n_epochs, 
    mean_toa_uncertainty=toa_uncertainty, 
    rms_residual=rms_residual, 
    observatory_code='@', 
    jitter_epoch=1, 
    output_file="simulated_toas_{0}.tim".format(n_epochs)
)
