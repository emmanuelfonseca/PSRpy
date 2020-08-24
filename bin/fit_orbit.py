#! /usr/bin/python

from PSRpy.orbit.velocities import doppler_shift_period
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import sys

plt.style.use('classic')

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

            a = input('Model parameters (ps, x, pb, ecc, om, t0): ')
            ps, x, pb, ecc, om, t0 = np.array(a.split(), dtype=float)
            mjds_model = np.linspace(np.min(mjds), np.max(mjds), num=1000)

            plt.errorbar(mjds, period, yerr=period_err, fmt='ro')
            initial_model = doppler_shift_period(ps, mjds_model, x, pb, ecc, om, t0)
            plt.plot(mjds_model, initial_model, 'b-')
            plt.grid()
            plt.show()

            repeat = input('Re-do model? (y/n): ')

            if (repeat == 'n'):
                keep_looping = False
                binary_pars[:] = ps, x, pb, ecc, om, t0
            else:
                pass

    return binary_pars


# read/unload data and create initial-guess model.
infile = (sys.argv)[1]
data = np.loadtxt(infile)
mjds = data[:, 0]
periods = data[:, 1]
periods_err = data[:, 2]
pars = guess_binary_model(mjds, periods, periods_err)

sys.exit()

# now, fit the data.
def chisq(parameters):
    """
    Computes the sum of square of differences between data and model of the data, 
    specifically for the fit_orbit.py routine in modeling spin-period variations.
    """

    p0, p1, p2, p3, p4, p5 = parameters
    model = p0 * (1 + doppler_shift_orbit(mjds, p1, p2, p3, p4, p5))
    return np.sum(((periods - model) / periods_err)**2) 

keep_looping = True

while keep_looping:

    result = minimize(chisq, pars, method='Nelder-Mead')
    print("Fit success: {0}".format(result.success))

    pars = result.x
    plt.subplot(211)
    plt.errorbar(mjds, periods, yerr=periods_err, fmt='ro')
    initial_model = pars[0] * (1 + doppler_shift_orbit(mjds, pars[1], pars[2], pars[3], pars[4], pars[5]))
    plt.plot(mjds, initial_model, 'b-')
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
print("    * PS = {0:.6f} ms".format(pars[0]))
print("    * A1 = {0:.6f} lt-s".format(pars[1]))
print("    * PB = {0:.6f} days".format(pars[2]))
print("    * E  = {0:.6f}".format(pars[3]))
print("    * OM = {0:.6f} deg".format(pars[4]))
print("    * T0 = {0:.6f} MJD".format(pars[5]))
