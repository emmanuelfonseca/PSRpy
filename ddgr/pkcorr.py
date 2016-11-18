#! /usr/bin/python

import astropy.constants as const
from numpy import sin, cos, sqrt
import numpy as np

# define some constants.
v_light = const.c.value
G     = const.G.value
Msun  = const.M_sun.value
Tsun  = G * Msun / v_light**3
d2r   = np.pi / 180
pc_in = 3.08567758e16
R0_in = 8.34           # distance to Galactic center, kpc.
R0err_in = 0.16           # distance error, kpc
v0_in = 240.           # circular speed of SSB, km/s
v0err_in = 8.             # circular speed error, km/s
c = 2.9979e8

def doppler(d,derr,b,l,mu,muerr):
    """
    Calculate Doppler-bias contribution (i.e. D-dot/D) to orbital decay 
    measurement.
    """
    # header stuff.
    R0, R0err = R0_in, R0err_in 
    v0, v0err = v0_in, v0err_in 
    pc = pc_in
    R0new, R0errnew = np.array([R0,R0err])*pc*1000.   # convert to m
    v0new, v0errnew = np.array([v0,v0err])*1000.      # convert to m/s
    z      = d*sin(b)
    d     *= pc*1000.
    derr  *= pc*1000.
    mu    *= 1. / 1000. / 3600. * d2r / 86400. / 365.25
    muerr *= 1. / 1000. / 3600. * d2r / 86400. / 365.25
    beta   = d/R0new*cos(b)-cos(l)
    den    = sin(l)**2+beta**2
    # expected bias from Galactic potential, diff. acceleration 
    # and Shklovskii effect.
    c1, c2, c3, c4 = -1.08e-19, 1.25, 0.0324, 0.58
    galpot = c1*(c2*z/sqrt(z**2+c3)+c4*z)*sin(b)
    galrel = -cos(b)*(v0new**2)/R0new*(cos(l)+beta/den)/c
    shklov = (mu**2)*d/c
    # now compute uncertainty.
    dgalpotdd = c1*sin(b)*sin(b)*(c2/sqrt(z**2+c3)-c2*z**2/(z**2+c3)**(1.5)+c4)/pc/1000.
    dgalreldd = (v0new*cos(b)/R0new)**(2)/den/c*(1.-2.*(beta**2)/den)
    dshklovdd = (mu**2)/c
    dshklovdm = 2.*mu*d/c
    err  = sqrt(((dgalpotdd+dgalreldd+dshklovdd)*derr)**2+(dshklovdm*muerr)**2)
    return np.array([galpot,galrel,shklov]), err


def distGR(xpbd,xpbderr,pb,b,l,mu,muerr,nmc=500):
    """
    Determine the distance to pulsar-binary system, assuming GR is correct.
    This algorithm uses a Monte Carlo method to obtain a distribution of 
    distances. 
    """
    # header stuff.
    pc = pc_in
    c = v_light / 1000 / pc
    R0, R0err = R0_in, R0err_in
    v0, v0err = v0_in, v0err_in
    pi = np.pi
    dist = np.zeros(nmc)
    R0   = np.random.normal(R0,R0err,size=(nmc))
    v0   = np.random.normal(v0,v0err,size=(nmc))/pc
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
            if (np.fabs(d-db) < 1e-5):
                dist[i] = d
                break
    v0 = v0*pc
    xpbd = xpbd*1e12
    mu = mu*1000.*3600./d2r*86400.*365.25
    return dist, R0, v0, mu, xpbd
