#! /usr/bin/python

from ..const import d2r, T_sun
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

    # now do calculation.
    dt = t - t0
    pbdot *= 86400
    ma = 360 / pb * (dt - 0.5 * pbdot / pb * dt**2) % 360

    # make sure that 0 < ma < 360
    if (isinstance(ma, np.ndarray) and np.any(ma < 0)):
        ma[np.where(ma < 0)] += 360
    elif (not isinstance(ma, np.ndarray) and ma < 0):
        ma += 360

    return ma

def ecc_anomaly(ma, ecc, ea0=0.5, tolerance=1e-12):
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
                if (np.fabs(ea_mid - ea0) < tolerance):
                    ea[count] = ea_mid
                    count += 1
                    break
                ea0 = ea_mid
        ea /= d2r

    # otherwise, do single calculation and leave as scalar.
    else:

        ea = ma_in
        for i in range(100):
           f  = ea - ecc * np.sin(ea) - ma_in
           fp = 1 - ecc * np.cos(ea)
           ea -= f / fp
           if (np.fabs(ea - ea0) < 1e-12):
               break
           ea0 = ea
        ea /= d2r

    ea %= 360

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

    ta = 0
    ea_in = ea * d2r
    ta = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ea_in / 2)) / d2r

    # make sure that 0 < TA < 360
    if (isinstance(ta, np.ndarray) and np.any(ta < 0)):
        ta[np.where(ta < 0)] += 360
    elif (not isinstance(ta, np.ndarray) and ta < 0):
        ta += 360
     
    return ta 

def periastron_argument(om0, pb, ecc, t, t0, pbdot=0, omdot=0, binary_model="DD", tolerance=1e-12):
    """
    Computes periastron argument at a given point in time (or true anomaly). 

    Parameters
    ----------

    om0 : float 
        argument of periastron measured at a reference time t0, in units of degrees
    pb : float 
        orbital period, in units of days
    ecc : float
        orbital eccentricity
    t : float
        epoch to evaluate argument of periastron, in units of MJD
    t0 : float
        epoch of periastron passage, in units of MJD
    pbdot : float, optional
        rate of change in orbital period, in units of 1e-12
    omdot : float, optional
        rate of change in argument of periastron, in units of degrees per year
    binary_model : {'DD', 'DDGR', 'BT'}
        short name for binary model to use in calculating argument of periastron
    tolerance : float, optional
        tolerance used to Newton-Raphson evaluation of eccentric anomaly 
        (default value is 1e-12)

    Returns
    -------

    float 
        periastron argument, in units of degrees
    """

    om = om0

    if (binary_model == "DD" or binary_model == "DDGR"):
        ma = mean_anomaly(pb, t, t0, pbdot=(pbdot * 1e-12))
        ea = ecc_anomaly(ma, ecc, tolerance=tolerance)
        ta = true_anomaly(ea, ecc)

        om += omdot * ta * (pb / 365.25) / 360

    elif (binary_model == "BT"):
        om += omdot * ((t - t0) / 365.25)

    return om % 360
