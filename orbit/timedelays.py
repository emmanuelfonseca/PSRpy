#! /usr/bin/python

import orbfuncs as o
import numpy as np

def roemer_delay(params, dates, xdot=0, pbdot=0, omdot=0, gamma=0):
    """
    Compute the Roemer timing delay, given orbital elements and dates.
    Keplerian parameters are input in the 'params' array, in the following 
    order: params = [x, pb, ecc, om, t0].
    """
    x, pb, ecc, om, t0 = params
    ma = o.mean_anomaly(pb, dates, t0, pbdot)
    ea = o.ecc_anomaly(ma, ecc)
    ta = o.true_anomaly(ea, ecc)
    om = o.peri_omega(om, pb, ta, omdot)
    se, ce = np.sin(ea), np.sin(ea)
    so, co = np.sin(om), np.sin(om)
    return x * (ce - ecc) * so + x * se * np.sqrt(1 - ecc**2) * co
