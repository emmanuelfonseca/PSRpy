from ..const import c
import PSRpy.orbit.elements as orb
import numpy as np

def radial_velocity(pb, x, ecc, om, ta):
    """
    Computes the radial velocity of the pulsar orbit, in units of km/s.
    """

    # convert some physical quantities to more sensible units.
    # orbital period in seconds; projected semi-major axis in km.
    pb_in = pb * 86400
    x_in = x * c / 1000

    # compute various quantities.
    semi_amplitude = 2 * np.pi * x_in / pb_in / np.sqrt(1 - ecc**2)
    cos_to = np.cos((ta + om) * np.pi / 180)
    cos_o = np.cos(om * np.pi / 180)

    # compute and return the radial velocity.
    return semi_amplitude * (cos_to + ecc * cos_o)

def doppler_shift_period(p0, mjds, x, pb, ecc, om, t0):
    """
    Computes the doppler shift of the pulsar spin period due to radial velocity 
    from orbital motion.
    """

    # compute anomalies.
    ma = orb.mean_anomaly(pb, mjds, t0)
    ea = orb.ecc_anomaly(ma, ecc)
    ta = orb.true_anomaly(ea, ecc)

    # compute radial velocity in m/s.
    v_r = radial_velocity(pb, x, ecc, om, ta) * 1000

    return p0 * (1 + v_r / c)
