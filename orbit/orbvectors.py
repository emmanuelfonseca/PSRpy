#! /usr/bin/pythin

from PSRpy.const import au, c, d2r
import orbfuncs as of
import numpy as np
import sys

def radius_eccentric_orbit(x, pb, ecc, om, t0, t, incl, asc, ecl_b=None, ecl_l=None, basis=1):
    """
    Returns the Cartesian coordinates of the orbital radius vector for a pulsar 
    in a binary system. The basis is taken to be the ecliptic coordinate system.

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
        - t = input epoch [MJD]
        - incl = system inclination [deg]
        - asc = longitude of ascending node [deg]
        - ecl_b = ecliptic latitude [deg]
        - ecl_l = ecliptic longitude [deg]
        - basis = basis of coordinate system:
            * 1 = plane of sky
            * 2 = ecliptic coordinate system

    Output:
        - radius vector [a.u., a.u., a.u.]
    """

    incl_in = incl * d2r
    asc_in = asc * d2r

    so, co = np.sin(asc_in), np.cos(asc_in)
    si, ci = np.sin(incl_in), np.cos(incl_in)

    ma = of.mean_anomaly(pb, t, t0)
    ea = of.ecc_anomaly(ma, ecc)
    ta = of.true_anomaly(ea, ecc)
    tot = (om + ta) * d2r

    # compute length of vector.
    ap = x / si * c / au
    r_mag = ap * (1 - ecc**2) / (1 + ecc * np.cos(ta * d2r))

    st, ct = np.sin(tot), np.cos(tot)

    X = -ap * (so * ct + co * ci * st)
    Y = ap * (co * ct - so * ci * st)
    Z = ap * (-si * st)

    # determine which basis to use, output results. 
    if (basis == 1):
        return np.array([X, Y, Z])
    elif (basis == 2):
        if (ecl_b is None or ecl_l is None):
            sys.exit("Error: must set ecliptic coordinates!")
        else:
            sb, cb = np.sin(ecl_b * d2r), np.cos(ecl_b * d2r)
            sl, cl = np.sin(ecl_l * d2r), np.cos(ecl_l * d2r)
            Xecl = X * sl - Y * cl * sb - Z * cl * cb
            Yecl = -(X * cl + Y * sl * sb + Z * sl * cb)
            Zecl = Y * cb - Z * sb
            return np.array([Xecl, Yecl, Zecl])
