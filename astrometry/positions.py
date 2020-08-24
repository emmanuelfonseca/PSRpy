#! /usr/bin/python

from astropy.coordinates import CartesianRepresentation
from ..const import d2r
import astropy.units as u
import numpy as np

def pulsar_position_ecliptic(ecl_latitude, ecl_longitude):
    """
    Returns the components of a unit vector pointing in the direction of a pulsar,  
    given its ecliptic coordinates. 

    Parameters
    ----------
    ecl_latitude : float 
        ecliptic latitude of the pulsar, in degrees.

    ecl_longitude : float
        ecliptic longitude of the pulsar, in degrees.

    Returns
    -------
    position_unit_vector : array_like
        unit vector pointing in the direction of the pulsar.
    """

    b = ecl_latitude * d2r
    l = ecl_longitude * d2r
    unit_x = np.cos(b) * np.cos(l)
    unit_y = np.cos(b) * np.sin(l)
    unit_z = np.sin(b)
    unit_xyz = [unit_x, unit_y, unit_z]

    # turn the vector into an Astropy-compatible object.
    position_unit_vector = CartesianRepresentation(unit_xyz, unit=u.dimensionless_unscaled)

    return position_unit_vector
