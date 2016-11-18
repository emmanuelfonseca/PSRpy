#! /usr/bin/python 

from const import d2r, spy, G, c, Msun, au
from numpy import sin, cos

def m2sini(fdots,m1,ecc,omg,truea,truead):
    """
    Compute expected mass of companion to within a factor of sine(inclination) 
    using the Josh & Rasio (1997) method for orbital determination from time 
    derivatives of the spin frequencies. See Equation 17 of JR97 and 
    'getorb' function.

    Required arguments:
        - 'fdots'  = array of first and second frequency derivatives.
        - 'm1'     = mass of pulsar (or inner-orbit objects, for a triple system; in Msun).
        - 'ecc'    = eccentricity.
        - 'omg'    = argument of periastron (in degrees).
        - 'truea'  = true anomaly (in degrees).
        - 'truead' = time-derivative of true anomaly (in degrees.)
    """
    # convert input to have necessary units.
    omg, truea = omg*d2r, truea*d2r
    truead /= spy
    # do the math, per Equation 17 in JR97.
    A   = 1. + ecc*cos(truea)
    m1 *= Msun 
    return -fdots[1]*c/fdots[0]/sin(truea+omg)*((m1**2)*(A**2)/G/truead**4)**(1./3.)/Msun



def a2(fdots,m1,ecc,omg,truea,truead):
    """
    Compute expected semimajor axis of companion using the Josh & Rasio (1997) method 
    for orbital determination from time derivatives of the spin frequencies. See 
    Equation 17 of JR97 and 'getorb' function.

    Required arguments:
        - 'fdots'  = array of first and second frequency derivatives.
        - 'm1'     = mass of pulsar (or inner-orbit objects, for a triple system; in Msun).
        - 'ecc'    = eccentricity.
        - 'omg'    = argument of periastron (in degrees).
        - 'truea'  = true anomaly (in degrees).
        - 'truead' = time-derivative of true anomaly (in degrees.)
    """
    m2s = m2sini(fdots,m1,ecc,omg,truea,truead)
    print m2s
    # exclude data points that correspond to eccentricity = 1.
    ecc    = ecc[0:len(ecc)-1]
    omg    = omg[0:len(omg)-1]
    truea  = truea[0:len(truea)-1]
    truead = truead[0:len(truead)-1]
    # convert input to have necessary units.
    omg, truea = omg*d2r, truea*d2r
    truead /= spy
    # do the math, per Equation 17 in JR97.
    A   = 1. + ecc*cos(truea)
    m1 *= Msun
    return -fdots[1]*m1*c*(A**2)/m2s/(1.-ecc**2)/fdots[0]/sin(truea+omg)/(truead**2)
