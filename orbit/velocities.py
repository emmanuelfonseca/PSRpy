from ..const import c
import numpy as np

def radial_velocity(pb, x, ecc, om, ta):
    """
    Computes the radial velocity of the pulsar orbit.
    """

    # convert some physical quantities to more sensible units.
    # orbital period in seconds; projected semi-major axis in km.
    pb_in = pb * 86400
    x_in = x * c / 1000

    # compute various quantities.
    semi_amplitude = 2 * np.pi * x / pb / np.sqrt(1 - ecc**2)
    cos_to = np.cos((ta + om) * np.pi / 180)
    cos_o = np.cos(om * np.pi / 180)

    # compute and return the radial velocity.
    return semi_amplitude * (cos_to + ecc * cos_o)
