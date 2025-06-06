#! /usr/bin/python

from PSRpy.const import c, pc, R0, R0_err, v0, v0_err
from numpy import pi, sin, cos, sqrt
import astropy.units as u
import numpy as np

d2r = pi / 180

def doppler(distance, distance_err, gal_b, gal_l, mu, mu_err, R0=R0.value, 
            R0_err=R0_err.value, v0=v0.value, v0_err=v0_err.value):
    """
    Calculates Doppler-bias contribution (i.e. D-dot/D) to orbital decay 
    measurement.

    Inputs
    ------

    distance : float
        Distance to pulsar, in units of kpc.

    distance_err : float
        Uncertainty in distance, in units of kpc.

    gal_b : float
        Galactic latitude, in units of radians.

    gal_l : float
        Galactic longitude, in units of radians.

    mu : float
        Composite proper motion, in units of milliarcseconds per year.

    mu_err : float
        Uncertainty in proper motion, in units of milliarcseconds per year.
    """

    # header stuff.
    R0new, R0_errnew = np.array([R0, R0_err]) * pc * 1000.   # convert to m
    v0new, v0_errnew = np.array([v0, v0_err]) * 1000.      # convert to m/s
    d = distance * pc * 1000.
    derr = distance_err * pc * 1000.
    z = distance * sin(gal_b)
    mu_in = mu / 1000. / 3600. * d2r / 86400. / 365.25
    mu_in_err = mu_err / 1000. / 3600. * d2r / 86400. / 365.25
    beta = d / R0new * cos(gal_b) - cos(gal_l)
    den = sin(gal_l)**2 + beta**2 

    # expected bias from Galactic potential, diff. acceleration, and Shklovskii effect.
    c1, c2, c3, c4 = -1.08e-19, 1.25, 0.0324, 0.58
    galpot = c1 * (c2 * z / sqrt(z**2 + c3) + c4 * z) * sin(gal_b)
    galrel = -cos(gal_b) * (v0new**2) / R0new * (cos(gal_l) + beta / den) / c.value
    shklov = (mu_in**2) * d / c.value

    # now compute uncertainty.
    dgalpotdd = c1 * sin(gal_b) * sin(gal_b) * (c2 / sqrt(z**2 + c3) - c2 * z**2 / \
                (z**2 + c3)**(1.5) + c4) / pc / 1000.
    dgalreldd = (v0new * cos(gal_b) / R0new)**(2) / den / c.value * (1. - 2. * (beta**2) / den)
    dshklovdd = (mu_in**2) / c.value
    dshklovdm = 2 * mu_in * d / c.value
    err = np.array([np.fabs(dgalpotdd) * derr, np.fabs(dgalreldd) * derr, 
          sqrt((dshklovdd * derr)**2 + (dshklovdm * mu_in_err)**2)])

    return np.array([galpot,galrel,shklov]), err


def distGR(xpbd, xpbderr, pb, gal_b, gal_l, mu, muerr, nmc=500, tolerance=1e-12):
    """
    Determine the distance to pulsar-binary system, assuming GR is correct.
    This algorithm uses a Monte Carlo method to obtain a distribution of 
    distances. 
    """
    # header stuff.
    dist = []
    R0mc = np.random.normal(R0.value, R0_err.value, size=nmc) * u.kpc
    v0mc = np.random.normal(v0.to(u.kpc / u.s).value, v0_err.to(u.kpc / u.s).value, size=nmc) * u.kpc / u.s    
    xpbdmc = np.random.normal(xpbd, xpbderr, size=nmc) * 1e-12 * u.s / u.s
    mumc = (np.random.normal(mu, muerr, size=nmc) * u.mas / u.yr).to(u.rad / u.s)
    pb = (pb * u.d).to(u.s)

    # compute trig terms ahead of time for convenience.
    cl = cos(gal_l * u.deg)
    sl = sin(gal_l * u.deg)
    cb = cos(gal_b * u.deg)
    sb = sin(gal_b * u.deg)

    # define Kuijken & Gilmore (1989) coefficients.
    c1 = -1.08e-19 / u.s
    c2 = 1.25 * u.kpc
    c3 = 0.0324 * u.kpc * u.kpc
    c4 = 0.58 * u.dimensionless_unscaled

    # now do the fun stuff, using Newton-Raphson method.
    for i in range(nmc):
        d = 1. * u.kpc

        for j in range(100):
            db = d
            z = d * sb
            beta = d / R0mc[i] * cb - cl  
            galpot = c1 * (c2 * z / sqrt(z**2 + c3) + c4 * z) * sb / u.kpc 
            galrel = cb * v0mc[i]**2 / R0mc[i] * (cl + beta / (sl**2 + beta**2)) / c.to(u.kpc / u.s)
            shklov = mumc[i]**2 * d / c.to(u.kpc / u.s) / u.rad**2
            # take derivatives of contributions - check these again!
            galpotd = c1 * (c2 * sb/sqrt(z**2 + c3) - \
                      c2 * z**2 * sb / (z**2 + c3)**(1.5) + c4) / u.kpc * sb
            galreld = cb * v0mc[i]**2 / c.to(u.kpc / u.s) / R0mc[i]**2 * (1 / (sl**2 + beta**2) - \
                      2 * beta**2 * cb / (sl**2 + beta**2)**2)
            shklovd = mumc[i]**2 / c.to(u.kpc / u.s) / u.rad**2
            f  = xpbdmc[i] + (galpot + galrel - shklov) * pb
            fp = (galpotd + galreld - shklovd) * pb
            d = d - f / fp
            if (np.fabs(d.value - db.value) < tolerance):
                dist.append(d.value)
                break

    return np.array(dist), R0mc.value, v0mc.value, mumc.to(u.mas / u.yr).value, (xpbdmc * 1e12).value
