#! /usr/bin/python

from .ssbfuncs import planet_position_ecliptic, pulsar_position_ecliptic
from .elements import anomaly_mean, anomaly_eccentric, anomaly_true, argument_periastron
from .orbvectors import radius_orbit_eccentric
from ..const import au, c, d2r, pc, T_sun
import numpy as np
import sys

def delay_orbit_parallax(time: float, period: float, axis_semimajor_projected: float, 
    eccentricity: float, argument_periastron: float, t0: float, inclination: float, 
    longitude_node_ascending: float, latitude_ecliptic: float, longitude_ecliptic: float):
    """
    Computes the parallax timing delay due to the orbital motion of the pulsar, given a set 
    of ecliptic coordinates, orbital parameters and distance measure.

    Parameters
    ----------
    time : array_like
        epoch at which to evaulate time delay 

    period : float 
        period of Keplerian orbit, in units of days

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    eccentricity : float
        Keplerian eccentricity for closed orbits (i.e., 0 <= eccentricity < 1)

    argument_periastron : float 
        Keplerian argument of periastron, in units of degrees

    t0 : float 
        epoch of periastron passage, in units of MJD

    inclination : float 
        geometric inclination of orbital system, in units of degrees

    longitude_node_ascending : float 
        geometric longitude of ascending node for orbit, in units of degrees

    latitude_ecliptic : float
        ecliptic longitude of binary system

    longitude_ecliptic : float
        ecliptic longitude of binary system

    Returns
    -------
    delay : array_like
        time delay due to orbital parallax motion
    """

    # compute terms needed to define radius vector.
    sini = np.sin(inclination * np.pi / 180)
    ap = axis_semimajor_projected / sini * c / au

    # now compute position and direction vectors.
    direction = ssb.pulsar_position_ecliptic(latitude_ecliptic, longitude_ecliptic)
    position_pulsar = ov.radius_orbit_eccentric(time, ap, orbital_period, eccentricity, 
        argument_periastron, t0, inclination, longitude_node_ascending, 
        latitude_ecliptic=latitude_ecliptic, longitude_ecliptic=longitude_ecliptic, ecliptic=True)
 
    # now compute cross product.
    delay = None

    if isinstance(epoch, np.ndarray):
        delay = np.zeros(len(time))

        for ii in range(len(time)):
            delay[ii] = np.sum(np.cross(position_pulsar[:, ii], direction)**2)

        delay = delay / 2 / c / (d * 1000 * pc)

    else:
        delay = np.sum(np.cross(position_pulsar, direction)**2) / 2 / c / (d * 1000 * pc)

    # now, return.
    return delay

def delay_orbit_annual_parallax(x, pb, ecc, om, t0, epoch, incl, asc, ecl_b, ecl_l, d=1, basis=2):
    """ 
    Computes the Roemer timing delay for the Solar System, given a set 
    of ecliptic coordinates.

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD]
        - epoch = epoch where delay is evaluated [MJD]
        - incl = system inclination [deg]
        - asc = longitude of ascending node [deg]
        - ecl_b = beta [deg]
        - ecl_l = lambda [deg]
        - d = distance [kpc]
        - basis = basis of coordinate system:
            * 1 = plane of sky
            * 2 = ecliptic coordinate system

    Output:
        - time delay [s]
    """

    # compute terms needed to define radius vector.
    sini = np.sin(inclination * np.pi / 180)
    ap = axis_semimajor_projected / sini * c / au

    # now define the position and direction vectors relevant for annual orbital parallax.
    direction = ssb.pulsar_position_ecliptic(latitude_ecliptic, longitude_ecliptic)
    position_earth = ssb.planet_position_ecliptic(epoch) * au
    position_pulsar = ov.radius_orbit_eccentric(time, ap, orbital_period, eccentricity, 
        argument_periastron, t0, inclination, longitude_node_ascending, 
        latitude_ecliptic=latitude_ecliptic, longitude_ecliptic=longitude_ecliptic, ecliptic=True)
 
    # now compute cross product.
    delay = None

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(time))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(np.cross(r_pulsar[:, ii], s) * np.cross(r_earth[:, ii], s))

        delay = delay / 2 / c / (d * 1000 * pc) 

    else:
        delay = np.sum(np.cross(r_earth, s) * np.cross(r_pulsar, s)) / 2 / c / (d * 1000 * pc)

    # now, return.
    return delay

def delay_orbit_correction_BT(dates, orbital_elements, xdot=0, pbdot=0, omdot=0, gamma=0, 
    eps1dot=0, eps2dot=0, m1=0, m2=0, dtheta=0, tolerance=1e-12, binary_model="DD", 
    orbital_phase=False):

    # first, compute Keplerian term.
    x, pb, ecc, om0, t0 = orbital_elements
    ma = anomaly_mean(dates, pb, t0, pbdot=(pbdot * 1e-12))
    ea = anomaly_eccentric(ma, ecc, nr_tolerance=tolerance)
    om = argument_periastron(dates, om0, pb, ecc, t0, pbdot=pbdot, omdot=omdot,
                               binary_model=binary_model, nr_tolerance=tolerance)
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    alpha = (x + (xdot * 1e-12) * (dates - t0)) * so
    beta = np.sqrt(1 - ecc**2) * (x + (xdot * 1e-12) * (dates - t0)) * co
    delay = alpha * (ce - ecc) + (beta + gamma) * se

    # next, compute second-order correction term.
    delay_corr = (alpha * se - beta * ce) * delay / ((pb * 86400) / 2 / np.pi) / (1 - ecc * ce)

    return delay_corr

def delay_orbit_roemer(dates, orbital_elements, xdot=0, pbdot=0, omdot=0, gamma=0, 
    eps1dot=0, eps2dot=0, m1=0, m2=0, dtheta=0, tolerance=1e-12, binary_model="DD", 
    orbital_phase=False):
    """
    Computes the Roemer timing delay for a pulsar in an eccentric binary system, given the 
    orbital elements and date(s). 

    Parameters
    ----------

    dates : array_like, float 
        epoch(s) to evaluate time delay, in units of MJD
    orbital_elements : tuple of floats
        the five orbital elements, in the following order for each listed binary model:
            - DD   : (x, pb, ecc, om, t0)
            - BT   : (x, pb, ecc, om, t0)
            - DDGR : (x, pb, ecc, om, t0)
            - ELL1 : (x, pb, eps1, eps2, tasc)
        a description of the elements is provded in the Notes section below
    xdot : float, optional
        rate of change in x, in units of 1e-12.
    pbdot : float, optional
        rate of change in pb, in units of 1e-12
    omdot : float, optional
        rate of change in om, in units of degree per year.
    gamma : float, optional 
        parameter for time dilation / gravitational redshift, in units of seconds.
    m1 : float, optional
        pulsar mass, in units of solar masses (used by 'DDGR', if set)
    m2 : float, optional
        companion mass, in units of solar masses (used by 'DDGR', if set)
    dtheta : float
        deformation of eccentricity in angular equation of motion
    tolerance : float, optional
        tolerance used to Newton-Raphson evaluation of eccentric anomaly 
        (default value is 1e-12)
    binary_model : {'DD', 'BT', 'DDGR','ELL1'}
        short name for binary model to use in calculating the Roemer delay

    Returns
    -------
    
    float
        Roemer time delay due to orbital motion, in seconds.

    Notes
    -----
    
    x   : projected semimajor axis of pulsar orbit, in units of light-seconds.
    pb  : orbital period, in units of days.
    ecc : orbital eccentricity, with 0 <= ecc < 1 for bound orbits.
    om  : argument of periastron, in units of degrees.
    t0  : epoch of periastron passage, in units of MJD

    The timing formula for each model can be found as follows:
        - BT: Blandford & Teukolsky (ApJ, 1976, 205, 580)
        - DD:
        - ELL1:
    """

    eccentric_orbit_models = ["BT", "DD", "DDGR", "Keplerian"]

    if (binary_model in eccentric_orbit_models):
        x, pb, ecc, om0, t0 = orbital_elements

        # first, compute Keplerian term.
        ma = anomaly_mean(dates, pb, t0, pbdot=(pbdot * 1e-12)) 
        ea = anomaly_eccentric(ma, ecc, nr_tolerance=tolerance) 
        om = argument_periastron(dates, om0, pb, ecc, t0, pbdot=pbdot, omdot=omdot, 
                                   binary_model=binary_model, nr_tolerance=tolerance) 
        se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
        so, co = np.sin(om * d2r), np.cos(om * d2r)
        alpha = (x + (xdot * 1e-12) * (dates - t0)) * so
        beta = np.sqrt(1 - ecc**2) * (x + (xdot * 1e-12) * (dates - t0)) * co
        delay = alpha * (ce - ecc) + (beta + gamma) * se

        # next, compute second-order correction term.

        if (binary_model == "BT"):
            delay += (alpha * se - beta * ce) * delay / ((pb * 86400) / 2 / np.pi) / (1 - ecc * ce)
        elif (binary_model == "Keplerian"):
            # for purely Keplerian motion, remove gamma term if supplied gamma is non-zero.
            delay -= gamma * se

        return delay

    elif (binary_model == "ELL1"):
        x, pb, eps1, eps2, tasc = orbital_elements

        # compute orbital elements and their variations from supplied derivatives.
        nb = 2 * np.pi / pb
        nbdot = -2 * np.pi / pb**2 * (pbdot * 1e-12)
        xnew = x + (xdot * 1e-12) * (dates - tasc)
        eps1new = eps1 + (eps1dot * 1e-12) * (dates - tasc)
        eps2new = eps2 + (eps2dot * 1e-12) * (dates - tasc)

        # determine OMDOT and EDOT from derivatives and values of the EPS1/EPS2 parameters.
        om = np.arctan2(eps1, eps2) * 180 / np.pi + 180
        so, co = np.sin(om * np.pi / 180), np.cos(om * np.pi / 180)
        ecc = np.sqrt(eps1**2 + eps2**2)
        a = np.array([[so, ecc * co], [co, -ecc * so]])
        b = np.array([eps1dot * 1e-12, eps2dot * 1e-12])
        edot, omdot = np.linalg.solve(a, b)

        # compute eccentricity terms for proper calculation of effective orbital frequency.
        t0 = tasc + (om * np.pi / 180) / (nb + (omdot * np.pi / 180 / 365.25 / 86400))
        nbbar = nb + (omdot * np.pi / 180 / 365.25 / 86400) - nbdot * (t0 - tasc)
        phi = (nbbar * (dates - tasc) + 0.5 * nbdot * (dates - tasc)**2) % (2 * np.pi)

        # compute delay.
        delay = xnew * (np.sin(phi) + 0.5 * eps2 * np.sin(2 * phi) - 0.5 * eps1 * np.cos(2 * phi))

        if orbital_phase:
            return (phi, delay)
        
        else:
            return delay

def delay_orbit_shapiro(dates, pb, ecc, om, t0, m2, sini, pbdot=0, omdot=0):
    """
    Computes the Shapiro timing delay for a pulsar-binary sustem, given the 
    orbital elements and dates. 

    Inputs:
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
        - m2 = companion mass [M_sun]
        - sini = sine of inclination angle [  ]
        - dates = epochs to evaluate delay [MJD]
        - pbdot = time derivative in pb [  ]
        - omdot = time derivative in om [deg / yr]

    Output:
        - time delay [s]
    """

    ma = anomaly_mean(dates, pb, t0, pbdot=pbdot) 
    ea = anomaly_eccentric(ma, ecc) 
    ta = anomaly_true(ea, ecc) 
    om = argument_periastron(dates, om, pb, ecc, t0, omdot=omdot) 
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    delay = -2 * T_sun * m2 * np.log(1 - ecc * ce - sini * ((ce - ecc) * so + se * np.sqrt(1 - ecc**2) * co))

    return delay
