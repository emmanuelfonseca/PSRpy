#! /usr/bin/python

from ..const import d2r, T_sun
import numpy as np
import sys

def anomaly_mean(time: float, orbital_period: float, t0: float, pbdot: float = 0.):
    """
    Computes the Keplerian mean anomaly, given orbital period and time parameters.

    Parameters
    ----------
    time : array_like
        desired time for evaluation of orbit, in units of MJD

    orbital_period : float
        period of the orbit, in units of MJD

    t0 : float
        epoch of passage through periastron, in units of MJD

    pbdot : float, optional 
        time derivative of orbital period, in units of second per second

    Returns
    -------
    anomaly_mean : array_like
    """

    # check if input time is a list or NumPy array.
    if isinstance(time, list):
        time = np.array(t)

    elif isinstance(t, np.ndarray):
        pass

    # now do calculation.
    dt = time - t0
    pbdot_in = pbdot * 86400.
    anomaly_mean = 360. / pb * (dt - 0.5 * pbdot_in / pb * dt**2)
    anomaly_mean = np.mod(anomaly_mean, 360.)

    return anomaly_mean

def anomaly_eccentric(anomaly_mean: float, eccentricity: float, initial_guess: float =0.5, 
    tolerance: float = 1e-12):
    """
    Computes the Keplerian eccentric anomaly, given mean anomaly and eccentricity. 
    This function uses a Newton-Raphson method to solve the transcendental equation
    that relates the eccentric and mean anomalies.

    Parameters
    ----------
    anomaly_mean : array_like
        Keplerian mean anomaly, in units of degrees

    eccentricity : float
        Keplerian eccentricity for closed orbits (i.e., 0 <= eccentricity < 1)

    initial_guess : float, optional 
        initial guess of eccentric anomaly for NR method.

    Returns
    -------
    anomaly_eccentric : array_like
        Keplerian eccentric anomaly, in units of degrees
    """

    ma_in = anomaly_mean * d2r    
    ea = None

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
                if (np.fabs(ea_mid - initial_guess) < tolerance):
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

    # ensure anomaly is within accepted bounds.
    ea = np.mod(ea, 360.)

    return ea 

def anomaly_true(anomaly_eccentric: float, eccentricity: float):
    """
    Computes true anomaly, given eccentric anomaly and eccentricity.

    Parameters
    ----------
    anomaly_eccentric : array_like
        Keplerian eccentric anomaly, in units of degrees.

    eccentricity : float
        Keplerian eccentricity for closed orbits (i.e., 0 <= eccentricity < 1)

    Returns
    -------
    anomaly_true : array_like
        Keplerian true anomaly, in units of degrees
    """

    ea_in = anomaly_eccentric * d2r
    term_ecc = np.sqrt((1 + eccentricity) / (1 - eccentricity))
    anomaly_eccentric = 2 * np.arctan(term_ecc * np.tan(ea_in / 2)) / d2r

    # ensure anomaly is within accepted bounds.
    anomaly_true = np.mod(anomaly_true, 360.)
     
    return anomaly_true

def argument_periastron(time: float, initial_periastron: float, orbital_period: float, 
    eccentricity: float, t0: float, pbdot: float = 0, omdot: float = 0, binary_model: str ="DD", 
    tolerance: float = 1e-12):
    """
    Computes periastron argument at a given point in time (or true anomaly). 

    Parameters
    ----------

    time : argument
        epoch to evaluate argument of periastron, in units of MJD

    initial_periastron : float 
        argument of periastron measured at a reference time t0, in units of degrees

    orbital_period : float 
        period of the orbit, in units of days

    eccentricity : float
        orbital eccentricity

    t0 : float
        epoch of periastron passage, in units of MJD

    pbdot : float, optional
        rate of change in orbital period, in units of second per second

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

    om = initial_periastron

    # use if/else to treat accounting of time derivatve in argument differently.
    if binary_model in ["DD", "DDGR"]:

        # first define anomalies.
        ma = anomaly_mean(time, orbital_period, t0, pbdot=pbdot)
        ea = anomaly_eccentric(ma, eccentricity, tolerance=tolerance)
        ta = anomaly_true(ea, eccentricity)

        # now compute variation term in Damour & Deruelle framework.
        om += omdot * ta * (pb / 365.25) / 360

    elif binary_model == "BT":

        # otherwise, just use Blandford & Teukolsky method, a first order Taylor expansion.
        om += omdot * ((time - t0) / 365.25)

    # ensure argument is within accepted bounds.
    om = np.mod(om, 360.)

    return om
