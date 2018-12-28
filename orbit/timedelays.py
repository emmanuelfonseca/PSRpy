#! /usr/bin/python

from PSRpy.const import au, c, d2r, pc, T_sun
import orbvectors as ov
import ssbfuncs as ssb
import orbfuncs as o
import mpmath as mp
import numpy as np

mp.mp.prec = 80

def roemer_delay_ssb(epoch, ecl_b, ecl_l):
    """
    Computes the Roemer timing delay about the Solar System Barycentre at the position the Earth, 
    given a set of ecliptic coordinates.

    Inputs
    ------

    epoch : float 
        Observation timestamp where delay is evaluated, in units of MJD.

    ecl_b : float 
        Ecliptic latitude, in units of degrees.

    ecl_l = : float
        Ecliptic longitude, in units of degrees.

    Output
    ------
    
        Time delay, in units of seconds.
    """

    r_earth = ssb.planet_position_ecliptic(epoch)
    s = ssb.pulsar_position_ecliptic(ecl_b, ecl_l)

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(epoch))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(r_earth[:, ii] * s)

        return delay * au / c

    else:
        return np.sum(r_earth * s) * au / c


def earth_parallax_delay(epoch, ecl_b, ecl_l, dist):
    """
    Computes the annual-parallax timing delay for the Earth, given a set 
    of ecliptic coordinates and distance measure.

    Inputs:
        - epoch = epoch where delay is evaluated [MJD]
        - ecl_b = beta [deg]
        - ecl_l = lambda [deg]
        - d = distance [kpc]

    Output:
        - time delay [s]
    """

    r_earth = ssb.planet_position_ecliptic(epoch) * au
    s = ssb.pulsar_position_ecliptic(ecl_b, ecl_l)

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(epoch))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(np.cross(r_earth[:, ii], s)**2)

        return delay * 1 / 2 / c / (dist * 1000 * pc)

    else:
        return np.sum(np.cross(r_earth, s)**2) * 1 / 2 / c / (dist * 1000 * pc)

def orbital_parallax_delay(x, pb, ecc, om, t0, epoch, incl, asc, ecl_b, ecl_l, d=1, basis=2):
    """
    Computes the parallax timing delay for the pulsar-binary, given a set 
    of ecliptic coordinates, orbital parameters and distance measure.

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

    r_pulsar = ov.radius_eccentric_orbit(x, pb, ecc, om, t0, epoch, incl, asc, 
        ecl_b, ecl_l, basis=basis)
    r_pulsar *= au
    s = ssb.pulsar_position_ecliptic(ecl_b, ecl_l)

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(epoch))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(np.cross(r_pulsar[:, ii], s)**2)

        return delay / 2 / c / (d * 1000 * pc)

    else:
        return np.sum(np.cross(r_pulsar, s)**2) / 2 / c / (d * 1000 * pc)

def annual_orbital_parallax_delay(x, pb, ecc, om, t0, epoch, incl, asc, ecl_b, ecl_l, d=1, basis=2):
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

    r_earth = ssb.planet_position_ecliptic(epoch) * au
    r_pulsar = ov.radius_eccentric_orbit(x, pb, ecc, om, t0, epoch, incl, asc, 
        ecl_b, ecl_l, basis=basis) * au
    s = ssb.pulsar_position_ecliptic(ecl_b, ecl_l)

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(epoch))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(np.cross(r_pulsar[:, ii], s) * np.cross(r_earth[:, ii], s))

        return delay / 2 / c / (d * 1000 * pc) 

    else:
        return np.sum(np.cross(r_earth, s) * np.cross(r_pulsar, s)) / 2 / c / (d * 1000 * pc)

def pulsar_roemer_delay_BT(dates, x, pb, ecc, om, t0, xdot=0, pbdot=0, omdot=0, gamma=0, 
    tolerance=1e-12, mp_version=False):
    """
    Computes the Roemer timing delay for a pulsar in an eccentric binary system, given the 
    orbital elements and date(s). The timing formula is taken from Blandford & Teukolsky 
    (ApJ, 1976, 205, 580).

    Parameters
    ----------

    dates : float 
        epoch(s) to evaluate time delay, in units of MJD
    x : float
        projected semimajor axis of pulsar orbit, in units of light-seconds.
    pb : float
        orbital period, in units of days.
    ecc : float 
        orbital eccentricity, with 0 <= ecc < 1 for bound orbits.
    om : float 
        argument of periastron, in units of degrees.
    t0 : float 
        epoch of periastron passage, in units of MJD
    xdot : float, optional
        rate of change in x, in units of 1e-12.
    pbdot : float, optional
        rate of change in pb, in units of 1e-12
    omdot : float, optional
        rate of change in periastron argument, in units of degree per year.
    gamma : float, optional 
        parameter for time dilation / gravitational redshift, in units of seconds.

    Returns
    -------
    
    float
        Roemer time delay due to orbital motion, in seconds.
    """

    if mp_version:
        
        if isinstance(dates, list):
            roemer_delay = []

            for ii in range(len(dates)):
                print ii
                ma = o.mean_anomaly(pb, dates[ii], t0, pbdot=pbdot, mp_version=mp_version) 
                ea = o.ecc_anomaly(ma, ecc, tolerance=tolerance, mp_version=mp_version) 
                ta = o.true_anomaly(ea, ecc, mp_version=mp_version) 
                om = o.peri_omega(om, pb, ta, omdot=omdot) 
                se, ce = mp.sin(ea * mp.pi / mp.mpf('180')), mp.cos(ea * mp.pi / mp.mpf('180'))
                so, co = mp.sin(om * mp.pi / mp.mpf('180')), mp.cos(om * mp.pi / mp.mpf('180'))
                roemer_delay.append(x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co)

            return roemer_delay

        else:

            ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot, mp_version=mp_version) 
            ea = o.ecc_anomaly(ma, ecc, tolerance=tolerance, mp_version=mp_version) 
            ta = o.true_anomaly(ea, ecc, mp_version=mp_version) 
            om = o.peri_omega(om, pb, ta, omdot=omdot, mp_version=mp_version) 
            se, ce = mp.sin(ea * mp.pi / mp.mpf('180')), mp.cos(ea * mp.pi / mp.mpf('180'))
            so, co = mp.sin(om * mp.pi / mp.mpf('180')), mp.cos(om * mp.pi / mp.mpf('180'))    

            return x * (ce - ecc) * so + x * se * mp.sqrt(mp.mpf('1') - mp.power(ecc, 2)) * co

    else:

        # first, compute Keplerian term.
        ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot) 
        ea = o.ecc_anomaly(ma, ecc, tolerance=tolerance) 
        ta = o.true_anomaly(ea, ecc) 
        om = o.peri_omega(om, pb, ta, omdot=omdot) 
        se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
        so, co = np.sin(om * d2r), np.cos(om * d2r)
        alpha = x * so
        beta = np.sqrt(1 - ecc**2) * x * co
        delay = alpha * (ce - ecc) + (beta + gamma) * se

        # next, compute second-order correction term.
        delay += (alpha * se - beta * ce) * delay / (pb / 2 / np.pi) / (1 - ecc * ce)

        return delay

def pulsar_roemer_delay_ELL1(dates, x, pb, eps1, eps2, tasc, xdot=0, pbdot=0, eps1dot=0, eps2dot=0,
    tolerance=1e-12):
    """

    """

def pulsar_shapiro_delay(dates, pb, ecc, om, t0, m2, sini, pbdot=0, omdot=0):
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

    ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot) 
    ea = o.ecc_anomaly(ma, ecc) 
    ta = o.true_anomaly(ea, ecc) 
    om = o.peri_omega(om, pb, ta, omdot=omdot) 
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    return -2 * T_sun * m2 * (1 - ecc * ce - sini * ((ce - ecc) * so + se * np.sqrt(1 - ecc**2) * co))
