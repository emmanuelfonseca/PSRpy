#! /usr/bin/python

from .ssb_data import ssb_planet_positions
from .elements import ecc_anomaly
from ..const import d2r
import numpy as np
import sys

def planet_position_ecliptic(epoch, planet='Earth'):
    """
    Returns the (X, Y, Z) coordinates of the selected planet, based on the 
    approximated orbital elements stored in the ssb_data subdirectory.

    Inputs:
        - epoch = time at which the position is to be evaluated [MJD]
        - planet = name of object to calculate position [str; default = Earth]

    Output:
        - position vector [a.u., a.u., a.u.]
    """

    Teph0 = 2451525.0 - 2400000.5 
    spp = ssb_planet_positions[planet]    

    sma = spp['sma'] + spp['sma-dot'] * (epoch - Teph0) / 36525
    ecc = spp['ecc'] + spp['ecc-dot'] * (epoch - Teph0) / 36525
    inc = spp['inc'] + spp['inc-dot'] * (epoch - Teph0) / 36525
    mlt = spp['mlt'] + spp['mlt-dot'] * (epoch - Teph0) / 36525
    per = spp['per'] + spp['per-dot'] * (epoch - Teph0) / 36525
    asc = spp['asc'] + spp['asc-dot'] * (epoch - Teph0) / 36525 

    # compute argument of perihelion, mean/eccentric anomalies.
    # all angular quanties are converted to radians hereafter.

    asc *= d2r
    inc *= d2r
    omega = (per - asc) * d2r
    mean_anom = (mlt - per) % 360
    #mean_anom *= d2r
    ecc_anom = ecc_anomaly(mean_anom, ecc) * d2r

    # compute planet's heliocentric coordinates, noting that z_helio = 0.

    x_helio = sma * (np.cos(ecc_anom) - ecc)
    y_helio = sma * np.sqrt(1 - ecc**2) * np.sin(ecc_anom)

    # finally, compute J2000 ecliptic coordinates, with x-axis aligned towards equinox.

    sa, ca = np.sin(asc), np.cos(asc)
    si, ci = np.sin(inc), np.cos(inc)
    so, co = np.sin(omega), np.cos(omega)

    x_ecl = (co * ca - so * sa * ci) * x_helio + (-so * ca - co * sa * ci) * y_helio
    y_ecl = (co * sa + so * ca * ci) * x_helio + (-so * sa + co * ca * ci) * y_helio
    z_ecl = (so * si) * x_helio + (co * si) * y_helio

    return np.array([x_ecl, y_ecl, z_ecl])

def pulsar_position_ecliptic(ecl_beta, ecl_lambda):
    """
    Returns the components of a unit vector pointing in the direction of a pulsar,  
    given its ecliptic coordinates. 

    Inputs:
        - ecl_beta   = beta [deg]
        - ecl_lambda = lambda [deg]

    Output:
        - Cartesian unit vector
    """

    b = ecl_beta * d2r
    l = ecl_lambda * d2r

    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])

