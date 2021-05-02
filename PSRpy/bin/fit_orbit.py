#! /usr/bin/python

from PSRpy.orbit.velocities import doppler_shift_period
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

# define some functions specific to this script.
def guess_binary_model(mjds, period, period_err):
    """
    Creates an interactive viewer of proposed model to use as a starting guess 
    for fitting the orbital elements from spin-period data.
    """

    binary_pars = np.zeros(6)
    keep_looping = True
    first = True

    while keep_looping:

        if first:

            plt.errorbar(mjds, period, yerr=period_err, fmt='ro')
            plt.grid()
            plt.show()

            first = False

        else:

            a = input('Model parameters (ps, pb, x, ecc, om, t0): ')
            ps, x, pb, ecc, om, t0 = np.array(a.split(), dtype=float)
            mjds_model = np.linspace(np.min(mjds), np.max(mjds), num=1000)

            plt.errorbar(mjds, period, yerr=period_err, fmt='ro')
            initial_model = doppler_shift_period(mjds_model, ps, pb, x, ecc, om, t0)
            plt.plot(mjds_model, initial_model, 'b-')
            plt.grid()
            plt.show()

            repeat = input('Re-do model? (y/n): ')

            if (repeat == 'n'):
                keep_looping = False
                binary_pars[:] = ps, pb, x, ecc, om, t0
            else:
                pass

    return binary_pars

def chisq(parameters, use_ELL1=False):
    """
    Computes the sum of square of differences between data and model of the data, 
    specifically for the fit_orbit.py routine in modeling spin-period variations.
    """

    p0, p1, p2, p3, p4, p5 = parameters

    # now compute model and residuals.
    model = doppler_shift_period(mjds, p0, p1, p2, p3, p4, p5, use_ELL1=use_ELL1)
    residuals = ((periods - model) / periods_err)**2
    chisq = np.sum(residuals)

    return chisq

# now define a simple argparse object.
parser = argparse.ArgumentParser(description="an interactive script for determining rough " + \
    "orbital parameters from a timeseries of pulsar-spin periods."
)

parser.add_argument(
    "input_file",
    type=str,
    help="an ASCII file containing period timeseries and uncertainties for each " + \
         "period; it is assumed that data are in whitespace-separated columns of: " + \
         "timestamp; period (in ms); and uncertainty on period (in ms)."
)

parser.add_argument(
    "--ELL1",
    action="store_true",
    default=False,
    dest="use_ELL1",
    help="If set, hold eccentricity and periastron argument fixed during fit for all " + \
         "other orbital parameters."
)

# TODO: add option to fit circular orbit (i.e., don't fit for eccentricity).

# grab data from command line.
args = parser.parse_args()
input_file = args.input_file
use_ELL1 = args.use_ELL1

# now define needed variables and determine initial guess for parameters.
data = np.loadtxt(input_file, usecols=(0, 1, 2))
mjds = data[:, 0]
periods = data[:, 1]
periods_err = data[:, 2]
pars = guess_binary_model(mjds, periods, periods_err)

# now, fit the data.
keep_looping = True

while keep_looping:

    result = minimize(chisq, pars, args=(use_ELL1), method='Nelder-Mead')
    print("Fit success: {0}".format(result.success))

    pars = result.x
    plt.subplot(211)
    plt.errorbar(mjds, periods, yerr=periods_err, fmt='ro')
    initial_model = doppler_shift_period(
        mjds, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]
    )

    # now plot results.
    mjds_interp = np.linspace(mjds.min(), mjds.max(), num=1000)
    model_interp = doppler_shift_period(
        mjds_interp, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]
    )

    plt.plot(mjds_interp, model_interp, 'b-')
    plt.ylabel('Spin Period (ms)')
    plt.grid()
    plt.subplot(212)
    plt.errorbar(mjds, periods - initial_model, yerr=periods_err, fmt='ro')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Residual (ms)')
    plt.grid()
    plt.show()

    repeat = input('Re-do fit? (y/n):')
    
    if (repeat == 'y'):
        pass
    else:
        keep_looping = False

print("Final parameter estimates:")
print("    * PS = {0:.15f} ms".format(pars[0]))
print("    * A1 = {0:.10f} lt-s".format(pars[1]))
print("    * PB = {0:.10f} days".format(pars[2]))
print("    * T0 = {0:.10f} MJD".format(pars[5]))

if use_ELL1:
    print("    * EPS1 = {0:.10f}".format(pars[3]))
    print("    * EPS2 = {0:.10f} deg".format(pars[4]))

else:
    print("    * E  = {0:.10f}".format(pars[3]))
    print("    * OM = {0:.10f} deg".format(pars[4]))
