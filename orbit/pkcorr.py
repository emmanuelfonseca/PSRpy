#! /usr/bin/python

from PSRpy.const import c, G, M_sun, pc, R0, R0_err, v0, v0_err
from numpy import sin, cos, sqrt
import numpy as np

d2r = np.pi / 180

def doppler(distance, distance_err, gal_b, gal_l, mu, mu_err):
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
    R0new, R0_errnew = np.array([R0, R0_err]) *pc * 1000.   # convert to m
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
    galrel = -cos(gal_b) * (v0new**2) / R0new * (cos(gal_l) + beta / den) / c
    shklov = (mu_in**2) * d / c

    # now compute uncertainty.
    dgalpotdd = c1 * sin(gal_b) * sin(gal_b) * (c2 / sqrt(z**2 + c3) - c2 * z**2 / \
                (z**2 + c3)**(1.5) + c4) / pc / 1000.
    dgalreldd = (v0new * cos(gal_b) / R0new)**(2) / den / c * (1. - 2. * (beta**2) / den)
    dshklovdd = (mu_in**2) / c
    dshklovdm = 2 * mu_in * d / c
    err = np.array([np.fabs(dgalpotdd) * derr, np.fabs(dgalreldd) * derr, 
          sqrt((dshklovdd * derr)**2 + (dshklovdm * mu_in_err)**2)])

    return np.array([galpot,galrel,shklov]), err


def distGR(xpbd,xpbderr,pb,b,l,mu,muerr,nmc=500):
    """
    Determine the distance to pulsar-binary system, assuming GR is correct.
    This algorithm uses a Monte Carlo method to obtain a distribution of 
    distances. 
    """
    # header stuff.
    c *= 1. / 1000. / pc
    pi = np.pi
    dist = np.zeros(nmc)
    R0   = np.random.normal(R0,R0_err,size=(nmc))
    v0   = np.random.normal(v0,v0_err,size=(nmc))/pc
    xpbd = np.random.normal(xpbd,xpbderr,size=(nmc))*1e-12
    mu   = np.random.normal(mu,muerr,size=(nmc))/1000./3600.*d2r/86400./365.25
    pb   = pb*86400.
    # now do the fun stuff, using Newton-Raphson method.
    for i in range(nmc):
        d = 1.
        for j in range(100):
            db = d
            z = d*sin(b)
            beta = d/R0[i]*cos(b)-cos(l)  
            galpot = 1.08e-19*(1.25*z/sqrt(z**2+0.0324)+0.58*z)*sin(b)
            galrel = cos(b)*(v0[i]**2)/R0[i]*(cos(l)+beta/(sin(l)**2+beta**2))/c
            shklov = (mu[i]**2)*d/c
            # take derivatives of contributions - check these again!
            galpotd = 1.08e-19*(1.25*sin(b)/sqrt(z**2+0.0324)-\
                      1.25*z**2*sin(b)/(z**2+0.0324)**(1.5)+0.58)*sin(b)
            galreld = (v0[i]**2)*cos(b)/c/R0[i]*(1./(sin(l)**2+beta**2)-\
                      2.*beta**2*cos(b)/R0[i]/(sin(l)**2+beta**2)**2)
            shklovd = (mu[i]**2)/c
            f  = xpbd[i]+(galpot+galrel-shklov)*pb
            fp = (galpotd+galreld-shklovd)*pb
            d = d-f/fp
            if (np.fabs(d-db) < 1e-12):
                dist[i] = d
                break
    v0 = v0*pc
    xpbd = xpbd*1e12
    mu = mu*1000.*3600./d2r*86400.*365.25
    return dist, R0, v0, mu, xpbd
