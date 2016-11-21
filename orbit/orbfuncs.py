#! /usr/bin/pythin

from PSRpy.const import d2r
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
    ma_in = ma * d2r
    ea = 0

    # use Newton-Raphson method to obtain best value of EA.
    if (isinstance(ma, np.ndarray)):
       count = 0
       ea = np.zeros(len(ma))

       # compute EA for each MA, separately.
       for ma0 in  ma_in:
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

    ea /= d2r

    # make sure that 0 < EA < 360
    if (isinstance(ea, np.ndarray) and np.any(ea < 0)):
        ea[np.where(ea < 0)] += 360
    elif (not isinstance(ea, np.ndarray) and ea < 0):
        ea += 360

    return ea 

def true_anomaly(ea, ecc):
    """
    Compute true anomaly, given eccentric anomaly and eccentricity.
    """
    ea_in = ea * d2r

    ta = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ea_in / 2)) / d2r

    # make sure that 0 < TA < 360
    if (isinstance(ta, np.ndarray) and np.any(ta < 0)):
        ta[np.where(ta < 0)] += 360
    elif (not isinstance(ta, np.ndarray) and ta < 0):
        ta += 360
     
    return ta 

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
    pb_in = pb / 365.25 # convert to years.
    om = omega0 + omdot * ta * pb / 2 / np.pi % 360
    return om
