#! /usr/bin/python

import astropy.constants as const
import numpy as np
import sys

# define some constants.
c    = const.c.value
G    = const.G.value
Msun = const.M_sun.value
Tsun = G*Msun/c**3
pi   = np.pi

def PBDOT_GR(m1,m2,pb,e):
    """
    Calculate orbital decay, as expected from GR.
    """
    # declare constant (i.e.non-mass) terms.
    pb = np.float(pb)*86400.
    e  = np.float(e)
    fe = (1.+(73./24.*e**2)+(37./96.*e**4))*(1.-e**2)**(-3.5)
    A  = -192.*pi/5.*(pb/2./pi)**(-5./3.)*fe*Tsun**(5./3.)
    # calculate!
    return A*m1*m2*(m1+m2)**(-1./3.)

def A1DOT_GR(m1, m2, pb, e):
    """
    Calculate rate of orbital decay for semi-major axis.
    Note: XDOT_GR = this rate x SINI.
    """
    # compute non-mass terms.
    fe = (1 + (73./24.) * e**2 + (37./96.) * e**4) / (1 - e**2)**(7./2.)
    A = -64. / 5. * (2 * np.pi * Tsun / pb)**2 * fe
    # calculate!
    return A * m1 * m2**2 / (m1 + m2)

def EDOT_GR(m1, m2, pb, e):
    """
    Calculate rate of circularization from GR.
    """
    nb = 2 * np.pi / pb
    eterm = (1 + (121./304.) * e**2) / (1 - e**2)**(5./2.)
    A = -304. / 15. * Tsun**(5./3.) * nb**(8./3.) * eterm
    return A * m1 * m2 / (m1 + m2)**(1./3.) * e

def DTHETA_GR(m1, m2, pb, e):
    """
    Calculate the DTHETA orbital-shape correction due to GR.
    """
    nb = 2 * np.pi / pb
    A = (Tsun * nb)**(2./3.)
    return A * (7. * m1**2 / 2. + 6 * m1 * m2 + 2 * m2**2) / (m1 + m2)**(4./3.)

def precession_GR(m1, m2, pb, e):
    """
    Calculate rate of geodetic precession of the spin axis.
    """
    A = (2 * np.pi / pb)**(5./3.) * Tsun**(2./3.) / (1 - e**2)
    return A * m2 * (4 * m1 + 3 * m2) / 2 / (m1 + m2)**(4./3.)

# compute terms due to rotational aberration.

def etaA(ps, pb, incl, e, eta, lamb):
    """
    Component of aberration-A term absorbed by other params. 
    (See Lorimer & Kramer, pg. 226., Equation 8.70)
    """
    si = np.sin(incl)
    se = np.sin(eta)
    sl = np.sin(lamb)
    return -ps / pb * se / sl / si / np.sqrt(1 - e**2)

