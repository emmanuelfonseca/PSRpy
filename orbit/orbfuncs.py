#! /usr/bin/pythin

import numpy as np
import sys

def mean_anomaly(pb, t, t0, pbdot=0):
    """
    Compute mean anomaly, given orbital pb and time parameters.
    """
    dt = t - t0
    pbdot *= 86400 * 365.25
    ma = 360 / pb * (dt - 0.5 * pbdot / pb * dt**2) % 360
    # make sure that 0 < ma < 360
    if (ma < 0):
        ma += 360
    return ma

def ecc_anomaly(ma, ecc, ea0=0.5):
    """
    Compute eccentric anomaly, given mean anomaly and eccentricity.
    """
    ea_new = ea0
    # use Newton-Raphson method to obtain best value of EA.
    for i in range(100):
        f  = ea_new - ecc * np.sin(ea_new) - ma
        fp = 1 - ecc * np.cos(ea_new)
        ea_new -= f / fp
        if (np.fabs(ea_new - ea0) < 1e-12):
            break
        ea0 = ea_new
    return ea_new * 180 / np.pi

def true_anomaly(ea, ecc):
    """
    Compute true anomaly, given eccentric anomaly and eccentricity.
    """
    ta = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ea / 2))
    return ta * 180 / np.pi

def peri_omega(omega0, pb, ta, omdot=0):
    """
    Compute periastron argument, given initial value, orbital pb, 
    true anomaly (ta) and periastron advance.

    Units of input:
        - [omega0] = degrees,
        - [pb] = days,
        - [ta] = degrees,
        - [omdot] = deg / yr
    """
    pb /= 365.25 # convert to years.
    return omega0 + omdot * ta * pb / 2 / np.pi

def roemer_delay(params, dates, xdot=0, pbdot=0, omdot=0, gamma=0):
    """
    Compute the Roemer timing delay, given orbital elements and dates.
    Keplerian parameters are input in the 'params' array, in the following 
    order: params = [x, pb, ecc, om, t0].
    """
    x, pb, ecc, om, t0 = params
    ma = mean_anomaly(pb, dates, t0, pbdot)
    ea = ecc_anomaly(ma, ecc)
    ta = true_anomaly(ea, ecc)
    om = peri_omega(om, pb, ta, omdot)
    se, ce = np.sin(ea), np.sin(ea)
    so, co = np.sin(om), np.sin(om)
    return x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co 
