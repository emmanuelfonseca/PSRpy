#! /usr/bin/python

from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, CartesianRepresentation
from astropy.time import Time
from ..const import au, c, d2r, pc
from .positions import pulsar_position_ecliptic
import astropy.units as u
import numpy as np
import sys

def roemer_delay(epoch, ecl_latitude, ecl_longitude, ephemeris='de432s'):
    """
    Computes the Roemer timing delay about the Solar System Barycentre at the 
    position the Earth, given a set of ecliptic coordinates for a given pulsar.

    Parameters
    ----------
    epoch : array_like, float 
        observation timestamp(s) where delay is evaluated, in MJD format.

    ecl_longitudeatitude : float 
        ecliptic latitude, in units of degrees.

    ecl_longitude = : float
        ecliptic longitude, in units of degrees.

    Returns
    -------
    delay: array_like, float
        time delay due to orbital motion of the Earth, in seconds.
    """

    # before computing anythin, set the ephemeris context to be value 
    # supplied at the function call.
    solar_system_ephemeris.set(ephemeris)
    time = Time(epoch, format='mjd')

    # now comoute the required position vectors.
    r_earth = get_body_barycentric('earth', time)
    s_pulsar = pulsar_position_ecliptic(ecl_latitude, ecl_longitude)

    # finally, compute and return the delay.
    delay = r_earth.dot(s_pulsar).to(u.m) / c

    return delay

def parallax_delay(epoch, ecl_latitude, ecl_longitude, distance, ephemeris='de432s'):
    """
    Computes the annual-parallax timing delay for the Earth, given a set 
    of ecliptic coordinates and distance measure.

    Parameters
    ----------
    epoch : array_like, float 
        observation timestamp(s) where delay is evaluated, in MJD format.

    ecl_latitude : float 
        ecliptic latitude, in units of degrees.

    ecl_longitude : float
        ecliptic longitude, in units of degrees.

    distance : float
        distance to the pulsar or pulsar-binary system, in units of kpc.

    Returns
    -------
    delay : array_like, dloat
        time delay due to parallax motion, in seconds.
    """

    # before computing anythin, set the ephemeris context to be value 
    # supplied at the function call.
    solar_system_ephemeris.set(ephemeris)
    time = Time(epoch, format='mjd')
    distance *= u.kpc

    # now comoute the required position vectors.
    r_earth = get_body_barycentric('earth', time)
    s_pulsar = pulsar_position_ecliptic(ecl_latitude, ecl_longitude) 
    r_cross_s = r_earth.cross(s_pulsar)

    # finally, compute and return the delay.
    delay = r_cross_s.dot(r_cross_s).to(u.m**2) / 2 / c / distance.to(u.m)

    return delay
