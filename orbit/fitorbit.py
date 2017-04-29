from scipy.optimize import minimize
from PSRpy.const import d2r
import matplotlib.pyplot as plt
import PSRpy.orbit as orb
import numpy as np
import sys

def doppler_shift_orbit(dates, x, pb, ecc, om, t0):
    """
    Computes fraction of period induced by the Dopper effect. 
    (See Equations 8.24, 8.25 of Lorimer & Kramer, 2005.)

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
        - dates  = epochs to evaluate delay [MJD]

    Output:
        - Doppler factor, v_orb/c [  ]
    """

    ma = orb.mean_anomaly(pb, dates, t0)
    ea = orb.ecc_anomaly(ma, ecc)
    ta = orb.true_anomaly(ea, ecc)
    om = orb.peri_omega(om, pb, ta)
    so, co = np.sin(om * d2r), np.cos(om * d2r)

    amp = 2 * np.pi * x / pb / 86400 / np.sqrt(1 - ecc**2)
    return amp * (np.cos((om + ta) * d2r) + ecc * co)

def fit_orbit(pars, mjds, periods, periods_errs, circular=False):
    """
    Applies the Doppler-shifted spin period model to a a set of period measurements and 
    fits for the orbital elements, assuming an eccentric orbit.

    Inputs:
        - pars = a list of the parameters to be fit:
            * pars[0] = spin period [ms]
            * pars[1] = projected semimajor axis [lt-s]
            * pars[2] = orbital period [days]
            * pars[3] = eccentricity [  ]
            * pars[4] = argument of periastron [deg]
            * pars[5] = epoch of periastron passage [MJD]
        - mjds = array of measurement epochs [MJD]
        - periods = measurements of pulsar-spin period [ms]
        - periods_errs = measurement uncertainties [ms]

    Outputs:
        - list of best-fit parameters, same order and units as input parameter list.
    """

    if (circular and len(pars) != 4):
        sys.exit("Circular orbit can only have 4 fit parameters, but {0} are given.".format(len(pars)))

    def chisq(parameters):

        if circular:
            p0, p1, p2, p3 = parameters
            model = p0 * (1 + doppler_shift_orbit(mjds, p1, p2, 0., 0., p3))
            return np.sum(((periods - model) / periods_errs)**2)

        else:

            p0, p1, p2, p3, p4, p5 = parameters
            model = p0 * (1 + doppler_shift_orbit(mjds, p1, p2, p3, p4, p5))
            return np.sum(((periods - model) / periods_errs)**2)

    result = minimize(chisq, pars, method='Nelder-Mead')
    print "Fit success: {0}".format(result.success)

    pars = result.x
    plt.subplot(211)
    plt.errorbar(mjds, periods, yerr=periods_errs, fmt='ro')
    if circular:
        mjds_mod = np.linspace(min(mjds), max(mjds), num=500)
        initial_model_mod = pars[0] * (1 + doppler_shift_orbit(mjds_mod, pars[1], pars[2], 0., 0., pars[3]))
        plt.plot(mjds_mod, initial_model_mod, 'b-')
    else:
        mjds_mod = np.linspace(min(mjds), max(mjds), num=500)
        initial_model_mod = pars[0] * (1 + doppler_shift_orbit(mjds_mod, pars[1], pars[2], pars[3], pars[4], pars[5]))
        plt.plot(mjds_mod, initial_model_mod, 'b-')
    plt.ylabel('Pulse Period (ms)')
    plt.grid()
    plt.subplot(212)
    initial_model = 0.
    if circular:
        initial_model = pars[0] * (1 + doppler_shift_orbit(mjds, pars[1], pars[2], 0., 0., pars[3]))
    else:
        initial_model = pars[0] * (1 + doppler_shift_orbit(mjds, pars[1], pars[2], pars[3], pars[4], pars[5]))
    plt.errorbar(mjds, periods - initial_model, yerr=periods_errs, fmt='ro')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Residual (ms)')
    plt.grid()
    plt.show()

    print "Final parameter estimates:"
    print "    * PS = {0:.9f} ms".format(pars[0])
    print "    * A1 = {0:.9f} lt-s".format(pars[1])
    print "    * PB = {0:.9f} days".format(pars[2])
    if circular:
        print "    * T0 = {0:.9f} MJD".format(pars[3])
    else:
        print "    * E  = {0:.9f}".format(pars[3])
        print "    * OM = {0:.9f} deg".format(pars[4])
        print "    * T0 = {0:.9f} MJD".format(pars[5])

    return pars
