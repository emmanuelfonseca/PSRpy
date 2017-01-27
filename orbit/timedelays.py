#! /usr/bin/python

from PSRpy.const import au, c, d2r
import ssbfuncs as ssb
import orbfuncs as o
import numpy as np

def roemer_delay_ssb(epoch, ecl_beta, ecl_lambda):
    """
    Computes the Roemer timing delay for the Solar System, given a set 
    of ecliptic coordinates.

    Inputs:
        - epoch      = epoch where delay is evaluated [MJD]
        - ecl_beta   = beta [deg]
        - ecl_lambda = lambda [deg]
    """

    r_earth = ssb.planet_position_ecliptic(epoch)
    s = ssb.pulsar_position_ecliptic(ecl_beta, ecl_lambda)

    if (isinstance(epoch, np.ndarray)):
        delay = np.zeros(len(epoch))

        for ii in range(len(epoch)):
            delay[ii] = np.sum(r_earth[:, ii] * s)

        return delay * au / c

    else:
        return np.sum(r_earth * s) * au / c


def roemer_delay(x, pb, ecc, om, t0, dates, xdot=0, pbdot=0, omdot=0, gamma=0):
    """
    Computes the Roemer timing delay for a pulsar-binary sustem, given the 
    orbital elements and dates. 

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
    """

    ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot) 
    ea = o.ecc_anomaly(ma, ecc) 
    ta = o.true_anomaly(ea, ecc) 
    om = o.peri_omega(om, pb, ta, omdot=omdot) 
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    return x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co
