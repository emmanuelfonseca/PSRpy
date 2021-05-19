from ..const import c, d2r
from . import elements as elem
import numpy as np
import sys

def velocity_radial(time: float, period: float, axis_semimajor_projected: float, 
    eccentricity: float, argument_periastron: float, t0: float):
    """
    Computes the radial velocity of the pulsar orbit, in units of km/s.

    Parameters
    ----------
    time : array_like
        desired time for evaluation of orbit, in units of MJD

    period : float
        period of the orbit, in units of days

    t0 : float
        epoch of passage through periastron, in units of MJD

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    eccentricity : float
        Keplerian eccentricity for closed orbits (i.e., 0 <= eccentricity < 1)

    argument_periastron : float 
        Keplerian argument of periastron, in units of degrees

    Returns
    -------
    velocity_radial : array_like
        the radial velocity of the pulsar's orbit
    """

    # compute anomalies for provided time and orbital elements.
    ma = elem.anomaly_mean(time, period, t0)
    ea = elem.anomaly_eccentric(ma, eccentricity)
    ta = elem.anomaly_true(ea, eccentricity)

    # convert some physical quantities to more sensible units.
    # orbital period in seconds; projected semi-major axis in km.
    pb_in = period * 86400
    x_in = axis_semimajor_projected * c.value / 1000

    # compute semi-amplitude and trig terms.
    semi_amplitude = 2 * np.pi * x_in / pb_in / np.sqrt(1 - eccentricity**2)
    ct = np.cos((ta + argument_periastron) * d2r)
    co = np.cos(argument_periastron * d2r)

    # compute and return the radial velocity.
    velocity_radial = semi_amplitude * (ct + eccentricity * co)

    return velocity_radial

def velocity_radial_ELL1(time: float, period: float, axis_semimajor_projected: float,
    epsilon1: float, epsilon2: float, t0: float):
    """
    Computes the radial velocity of the pulsar orbit in a nearly circular orbit, 
    in units of km/s. The form is taken from the ELL1 timing model or orbital 
    motion; see A6 in Lange et al. (2001, MNRAS, 326, 274).

    Parameters
    ----------
    time : array_like
        desired time for evaluation of orbit, in units of MJD

    period : float
        period of the orbit, in units of days

    t0 : float
        epoch of passage through periastron, in units of MJD

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    epsilon1 : float
        the first Laplace-Lagrange parameter of the ELL1 model

    epsilon2 : float
        the second Laplace-Lagrange parameter of the ELL1 model

    Returns
    -------
    velocity_radial : array_like
        the radial velocity of the pulsar's orbit
    """

    # compute trig terms.
    frequency_orbit = 2 * np.pi / period 
    phase_orbit = frequency_orbit * (time - t0)
    cp = np.cos(phase_orbit)
    s2p = np.sin(2 * phase_orbit)
    c2p = np.cos(2 * phase_orbit)

    # compute and return the radial velocity.
    amplitude = axis_semimajor_projected * frequency_orbit * c.value / 1000 / 86400
    velocity_radial = amplitude * (cp + epsilon1 * c2p + epsilon2 * s2p)

    return velocity_radial

def doppler_shift_period(time: float, period_spin: float, period:float, 
    axis_semimajor_projected: float, eccentricity1: float, eccentricity2: float, 
    t0: float, use_ELL1: bool = False):
    """
    Computes the doppler shift of the pulsar spin period due to radial velocity 
    from orbital motion.

    Parameters
    ----------
    time : array_like
        desired time for evaluation of orbit, in units of MJD

    period : float
        period of the orbit, in units of days

    t0 : float
        epoch of passage through periastron, in units of MJD

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    eccentricity1 : float
        one of two 'eccentricity parameters' that define the shape and orietnation 
        of the orbit. if use_ELL1 = False, then this is the first Laplace-Lagrange 
        parameter; otherwise, this is the Keplerian eccentricity.

    eccentricity2 : float
        one of two 'eccentricity parameters' that define the shape and orietnation 
        of the orbit. if use_ELL1 = False, then this is the second Laplace-Lagrange 
        parameter; otherwise, this is the Keplerian argument of periastron..

    Returns
    -------
    velocity_radial : array_like
        the radial velocity of the pulsar's orbit
    """

    # compute radial velocity in m/s.
    v_r = np.zeros(len(time))

    if use_ELL1:
        v_r = velocity_radial_ELL1(
            time, period, axis_semimajor_projected, eccentricity1, eccentricity2, t0
        )

    else:
        v_r = velocity_radial(
            time, period, axis_semimajor_projected, eccentricity1, eccentricity2, t0
        )

    v_r *= 1000

    # now compute doppler-shifted velocity.
    period_shifted = period_spin * (1 + v_r / c.value)

    return period_shifted 
