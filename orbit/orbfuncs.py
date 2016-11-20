#! /usr/bin/pythin

import numpy as np
import sys

def mean_anomaly(pb, t, t0, pbdot=0):
    """
    Compute mean anomaly, given orbital pb and time parameters.
    """
    # check if input time is a list or NumPy array.
    if (isinstance(t, list)):
        t = np.array(t)
    elif (isinstance(t, np.ndarray)):
        pass

    # now do calculation
    dt = t - t0
    pbdot *= 86400. 
    ma = 360 / pb * (dt - 0.5 * pbdot / pb * dt**2) % 360


    # make sure that 0 < ma < 360
    if (isinstance(ma, np.ndarray) and np.any(ma < 0)):
        ma[np.where(ma < 0)] += 360
    elif (not isinstance(ma, np.ndarray) and ma < 0):
        ma += 360

    return ma

def ecc_anomaly(ma, ecc, ea0=0.5):
    """
    Compute eccentric anomaly, given mean anomaly and eccentricity.
    """
    ma *= np.pi / 180
    ea = 0

    # use Newton-Raphson method to obtain best value of EA.
    if (isinstance(ma, np.ndarray)):
       count = 0
       ea = np.zeros(len(ma))

       # compute EA for each MA, separately.
       for ma0 in  ma:
           ea_mid = ma0
           for i in range(100):
               f  = ea_mid - ecc * np.sin(ea_mid) - ma0
               fp = 1 - ecc * np.cos(ea_mid)
               ea_mid -= f / fp
               if (np.fabs(ea_mid - ea0) < 1e-12):
                   ea[count] = ea_mid
                   count += 1
                   break
               ea0 = ea_mid
    else:
        ea = ma
        for i in range(100):
           f  = ea - ecc * np.sin(ea) - ma
           fp = 1 - ecc * np.cos(ea)
           ea -= f / fp
           if (np.fabs(ea - ea0) < 1e-12):
               break
           ea0 = ea

    # make sure that 0 < EA < 360
    if (isinstance(ea, np.ndarray) and np.any(ea < 0)):
        ea[np.where(ea < 0)] += 360
    elif (not isinstance(ea, np.ndarray) and ea < 0):
        ea += 360

    return ea * 180 / np.pi

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
