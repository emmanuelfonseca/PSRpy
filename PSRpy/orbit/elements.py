#! /usr/bin/python

from ..const import d2r, T_sun
import numpy as np
import sys

def anomaly_mean(time: float, period: float, t0: float, pbdot: float = 0.):
    """
    Computes the Keplerian mean anomaly, given orbital period and time parameters.

    Parameters
    ----------
    time : array_like
        desired time for evaluation of orbit, in units of MJD

    period : float
        period of the orbit, in units of days

    t0 : float
        epoch of passage through periastron, in units of MJD

    pbdot : float, optional 
        time derivative of orbital period, in units of second per second

    Returns
    -------
    anomaly_mean : array_like
    """

    # ensure time variable is in NumPy ndArray format.
    dt = np.array(time) - t0

    # now perform calculation.
    anomaly_mean = 360. / period * (dt - 0.5 * pbdot / period * dt**2)
    anomaly_mean = np.mod(anomaly_mean, 360.)

    return anomaly_mean

def anomaly_eccentric(anomaly_mean: float, eccentricity: float, initial_guess: float = 0.5, 
    nr_attempts: int = 100, nr_tolerance: float = 1e-12):
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

    nr_attempts : int, optional
        number of iterations for attempting the Newton-Raphson method 

    nr_tolerance : float, optional
        tolerance used for evaluating accuracy of the Newton-Raphson method

    Returns
    -------
    anomaly_eccentric : array_like
        Keplerian eccentric anomaly, in units of degrees
    """

    
    # ensure time variable is in NumPy ndArray format, convert to radians.
    ma_in = np.array(anomaly_mean) * d2r    

    # define variables needed for NR method.
    ea = 0.
    ea0 = initial_guess
    ea_mid = ma_in.copy()

    # loop over desired number of attempts to perform NR calculation.
    for idx in range(nr_attempts):
        func = ea_mid - eccentricity * np.sin(ea_mid) - ma_in
        func_deriv = 1 - eccentricity * np.cos(ea_mid)
        ea_mid -= func / func_deriv

        # if the anomaly has sub-threshold accuracy, break out.
        if np.all(np.fabs(ea_mid - ea0) < nr_tolerance):
            ea = ea_mid
            break

        # otherwise, update prior-step array and redo calculation.
        else:
            ea0 = ea_mid

    # ensure anomaly is within accepted bounds.
    ea /= d2r
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

    
    # ensure time variable is in NumPy ndArray format, convert to radians.
    ea_in = np.array(anomaly_eccentric) * d2r

    # now calculate true anomaly.
    term_ecc = np.sqrt((1 + eccentricity) / (1 - eccentricity))
    ta = 2 * np.arctan(term_ecc * np.tan(ea_in / 2)) / d2r

    # ensure anomaly is within accepted bounds.
    ta = np.mod(ta, 360.)
     
    return ta

def argument_periastron(time: float, initial_periastron: float, period: float, 
    eccentricity: float, t0: float, pbdot: float = 0, omdot: float = 0, 
    binary_model: str ="DD", nr_attempts: int = 100, nr_tolerance: float = 1e-12):
    """
    Computes periastron argument at a given point in time (or true anomaly). 

    Parameters
    ----------

    time : argument
        epoch to evaluate argument of periastron, in units of MJD

    initial_periastron : float 
        argument of periastron measured at a reference time t0, in units of degrees

    period : float 
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

    nr_attempts : int, optional
        number of iterations for attempting the Newton-Raphson method 

    nr_tolerance : float, optional
        tolerance used for evaluating accuracy of the Newton-Raphson method

    Returns
    -------

    float 
        periastron argument, in units of degrees
    """

    om = np.zeros(len(np.array(time))) + initial_periastron

    # use if/else to treat accounting of time derivatve in argument differently.
    if binary_model in ["DD", "DDGR"]:

        # first define anomalies.
        ma = anomaly_mean(time, period, t0, pbdot=pbdot)
        ea = anomaly_eccentric(ma, eccentricity, nr_attempts=nr_attempts, nr_tolerance=nr_tolerance)
        ta = anomaly_true(ea, eccentricity)

        # now compute variation term in Damour & Deruelle framework.
        om += omdot * ta * (period / 365.25) / 360

    elif binary_model == "BT":

        # otherwise, just use Blandford & Teukolsky method, a first order Taylor expansion.
        dt = np.array(time) - t0
        dt /= 365.25 
        om += omdot * dt

    # ensure argument is within accepted bounds.
    om = np.mod(om, 360.)

    return om
