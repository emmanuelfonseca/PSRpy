#! /usr/bin/python

from .elements import anomaly_mean, anomaly_eccentric, anomaly_true
from ..const import au, c, d2r
import numpy as np
import sys

def radius_orbit_eccentric(time: float, axis_semimajor: float, orbital_period: float, 
    eccentricity: float, argument_periastron: float, t0: float, inclination: float, 
    longitude_node_ascending: float, latitude_ecliptic: float = None, longitude_ecliptic: float = None, 
    ecliptic: bool = False):
    """
    Returns the Cartesian coordinates of the orbital radius vector for a pulsar 
    in an eccentric binary system. The basis is taken to be the ecliptic coordinate system.

    Parameters
    ----------
    time : array_like
        epoch to evaluate orbital radius vector, in units of MJD

    axis_semimajor_projected : float 
        Keplerian semimajor axis, in units of AU

    orbital_period : float 
        period of Keplerian orbit, in units of days

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

    latitude_ecliptic : float, optional
        ecliptic longitude of binary system

    longitude_ecliptic : float, optional
        ecliptic longitude of binary system

    ecliptic : bool, optional 
        if True, compute vector with basis in the ecliptic coordinate system

    Returns
    -------
    vector : array_like
        three-dimension position vector for orbit, in units of AU
    """

    incl_in = inclination * d2r
    asc_in = longitude_node_ascending * d2r

    so, co = np.sin(asc_in), np.cos(asc_in)
    si, ci = np.sin(incl_in), np.cos(incl_in)

    ma = anomaly_mean(time, orbital_period, t0)
    ea = anomaly_eccentric(ma, eccentricity)
    ta = anomaly_true(ea, eccentricity)
    tot = np.mod(argument_periastron + ta, 360.) * d2r

    # compute length of vector.
    r_mag = axis_semimajor * (1 - eccentricity**2) / (1 + eccentricity * np.cos(ta * d2r))

    st, ct = np.sin(tot), np.cos(tot)

    X = -ap * (so * ct + co * ci * st)
    Y = ap * (co * ct - so * ci * st)
    Z = ap * (-si * st)

    # determine which basis to use, output results. 
    if ecliptic and all([latitude_ecliptic, longitude_ecliptic]):
        sb, cb = np.sin(ecl_b * d2r), np.cos(ecl_b * d2r)
        sl, cl = np.sin(ecl_l * d2r), np.cos(ecl_l * d2r)
        Xecl = X * sl - Y * cl * sb - Z * cl * cb
        Yecl = -(X * cl + Y * sl * sb + Z * sl * cb)
        Zecl = Y * cb - Z * sb
        X = Xecl
        Y = Yecl
        Z = Zecl

    elif ecliptic:
        print("WARNING: ecliptic basis is desired but one or both ecliptic coordinates are not set!")
        print("... returning vector with basis relative to plane of sky ...")

    # now, return vector.
    vector = np.array([X, Y, Z])

    return vector
