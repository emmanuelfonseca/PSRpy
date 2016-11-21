#! /usr/bin/python

from PSRpy.const import d2r
import orbfuncs as o
import numpy as np

def roemer_delay(params, dates, xdot=0, pbdot=0, omdot=0, gamma=0):
    """
    Compute the Roemer timing delay, given orbital elements and dates.
    Keplerian parameters are input in the 'params' array, in the following 
    order: params = [x, pb, ecc, om, t0].
    """
    x, pb, ecc, om, t0 = params
    ma = o.mean_anomaly(pb, dates, t0, pbdot=pbdot) 
    ea = o.ecc_anomaly(ma, ecc) 
    ta = o.true_anomaly(ea, ecc) 
    om = o.peri_omega(om, pb, ta, omdot=omdot) 
    se, ce = np.sin(ea * d2r), np.cos(ea * d2r)
    so, co = np.sin(om * d2r), np.cos(om * d2r)
    return x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co
