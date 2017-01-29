#! /usr/bin/pythin

from PSRpy.const import d2r
import numpy as np
import sys

def mean_anomaly(pb, t, t0, pbdot=0):
    """
    Computes mean anomaly, given orbital period and time parameters.

    Inputs:
        - pb = orbital period [days]
        - t = epoch to evaluate angle [MJD]
        - t0 = reference epoch [MJD]
        - pbdot = time derivative of orbital period [  ]

    Output:
        - mean anomaly [deg]
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
    Computes eccentric anomaly, given mean anomaly and eccentricity.

    Inputs:
        - ma = mean anomaly [deg]
        - ecc = orbital eccentricity [  ]
        - ea0 = initial guess of eccentric anomaly [  ]

    Output:
        - eccentric anomaly [deg]
    """

    ma_in = ma * d2r
    ea = 0

    # if MA is an array, loop over each entry and apply N-R method.
    if (isinstance(ma, np.ndarray)):
       count = 0
       ea = np.zeros(len(ma))

       # in this case, turn 'ecc' into an array.
       if (not isinstance(ecc, np.ndarray)):
           ecc = np.zeros(len(ma)) + ecc

       # compute EA for each MA, separately.
       for ma0, ecc0 in zip(ma_in, ecc):
           ea_mid = ma0
           for i in range(100):
               f  = ea_mid - ecc0 * np.sin(ea_mid) - ma0
               fp = 1 - ecc0 * np.cos(ea_mid)
               ea_mid -= f / fp
               if (np.fabs(ea_mid - ea0) < 1e-12):
                   ea[count] = ea_mid
                   count += 1
                   break
               ea0 = ea_mid

    # otherwise, do single calculation and leave as scalar.
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
    Computes true anomaly, given eccentric anomaly and eccentricity.

    Inputs:
        - ea = eccentric anomaly [deg]
        - ecc = orbital eccentricity [deg]

    Output:
        - true anomaly [deg]
    """

    ea_in = ea * d2r
    ta = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ea_in / 2)) / d2r

    # make sure that 0 < TA < 360
    if (isinstance(ta, np.ndarray) and np.any(ta < 0)):
        ta[np.where(ta < 0)] += 360
    elif (not isinstance(ta, np.ndarray) and ta < 0):
        ta += 360
     
    return ta 

def peri_omega(om0, pb, ta, omdot=0):
    """
    Computes periastron argument as a function of time, given initial 
    value, orbital pb, true anomaly (ta) and periastron advance.

    Input:
        - om0 = periastron argument [deg]
        - pb = orbital period [days]
        - ta = true anomaly [deg]
        - omdot = time derivative of periastron argument [deg / yr]

    Output:
        - periastron argument [deg]
    """
    pb_in = pb / 365.25 # convert to years.
    om = om0 + omdot * ta * pb / 2 / np.pi % 360
    return om
