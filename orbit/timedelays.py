#! /usr/bin/python

from PSRpy.const import au, c, d2r, pc
import ssbfuncs as ssb
import orbvectors as ov
import orbfuncs as o
import numpy as np

def roemer_delay_ssb(epoch, ecl_b, ecl_l):
    """
    Computes the Roemer timing delay for the Earth, given a set 
    of ecliptic coordinates.

    Inputs:
        - epoch      = epoch where delay is evaluated [MJD]
        - ecl_b   = beta [deg]
        - ecl_l = lambda [deg]

    Output:
        - time delay [s]
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


def earth_parallax_delay(epoch, ecl_b, ecl_l, d=1):
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

        return delay / 2 / c / (d * 1000 * pc)

    else:
        return np.sum(np.cross(r_earth, s)**2) / 2 / c / (d * 1000 * pc) 

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
        return np.sum(np.cross(r_earth, s)**2) / 2 / c / (d * 1000 * pc)

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

def pulsar_roemer_delay(x, pb, ecc, om, t0, dates, xdot=0, pbdot=0, omdot=0, gamma=0):
    """
    Computes the Roemer timing delay for a pulsar-binary sustem, given the 
    orbital elements and dates. 

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
        - dates  = epochs to evaluate delay [MJD]
        - xdot = time derivatve in x [  ]
        - pbdot = time derivative in pb [  ]
        - omdot = time derivative in om [deg / yr]
        - gamma = parameter for time dilation / gravitational redshift [s]

    Output:
        - time delay [s]
    """

    ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot) 
    ea = o.ecc_anomaly(ma, ecc) 
    ta = o.true_anomaly(ea, ecc) 
    om = o.peri_omega(om, pb, ta, omdot=omdot) 
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    return x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co
