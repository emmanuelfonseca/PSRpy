#! /usr/bin/python

from PSRpy.const import c, G, M_sun, T_sun
import numpy as np

pi   = np.pi

def gamma_GR(m1, m2, pb, e):
    """
    Calculates time-averaged gravitational redshift / time dilation parameter, 
    as expected from GR.
    """
    
    pb_in = pb * 86400
    A = e * (pb_in / 2 / pi)**(1./3.) * T_sun**(2./3.)
    return A * m2 * (m1 + 2 * m2) / (m1 + m2)**(4./3.)

def omdot_GR(m1, m2, pb, e):
    """
    Calculates periastron advance, as expected from GR.
    """

    pb_in = pb * 86400
    omdot = 3 * (pb_in / 2 / pi)**(-5./3.) * (T_sun * (m1 + m2))**(2./3.) / (1 - e**2)
    return omdot * 180 / pi * 86400 * 365.25

def pbdot_GR(m1, m2, pb, e):
    """
    Calculates orbital decay of binary period, as expected from GR.
    """

    pb_in = np.float(pb)*86400.
    fe = (1 + (73. / 24. * e**2) + (37. / 96. * e**4)) * (1 - e**2)**(-3.5)
    A  = -192 * pi / 5 * (pb_in / 2 / pi)**(-5./3.) * fe * T_sun**(5./3.)
    return A*m1*m2*(m1+m2)**(-1./3.)

def r_GR(m2):
    """
    Calculates Shapiro 'range' parameter, as expected from GR.
    """

    return T_sun * m2

def s_GR(m1, m2, pb, x):
    """
    Calculates Shapiro 'shape' parameter, as expected from GR / mass function.
    """

    pb_in = pb * 86400
    A = x * (pb_in / 2 / pi)**(-2./3.) * T_sun**(-1./3.)
    return A * (m1 + m2)**(2./3.) / m2

def xdot_GR(m1, m2, pb, e):
    """
    Calculates rate of orbital decay for semi-major axis.
    Note: XDOT_GR = this rate x SINI.
    """

    pb_in = pb * 86400
    fe = (1 + (73./24.) * e**2 + (37./96.) * e**4) / (1 - e**2)**(7./2.)
    A = -64. / 5. * (2 * pi * T_sun / pb_in)**2 * fe
    return A * m1 * m2**2 / (m1 + m2)

def edot_GR(m1, m2, pb, e):
    """
    Calculates rate of circularization from GR.
    """
    nb = 2 * np.pi / pb
    eterm = (1 + (121./304.) * e**2) / (1 - e**2)**(5./2.)
    A = -304. / 15. * T_sun**(5./3.) * nb**(8./3.) * eterm
    return A * m1 * m2 / (m1 + m2)**(1./3.) * e

def dtheta_GR(m1, m2, pb, e):
    """
    Calculate the dtheta orbital-shape correction due to GR.
    """
    nb = 2 * np.pi / pb
    A = (T_sun * nb)**(2./3.)
    return A * (7. * m1**2 / 2. + 6 * m1 * m2 + 2 * m2**2) / (m1 + m2)**(4./3.)

def precession_GR(m1, m2, pb, e):
    """
    Calculate rate of geodetic precession of the spin axis.
    """
    A = (2 * np.pi / pb)**(5./3.) * T_sun**(2./3.) / (1 - e**2)
    return A * m2 * (4 * m1 + 3 * m2) / 2 / (m1 + m2)**(4./3.)

# compute orthometric terms.

def stig_GR(sini):
    """
    Calculates orthometric parameterization of inclination.
    """

    cosi = np.sqrt(1 - sini**2)
    return sini / (1 + cosi)

def h3_GR(m2, sini):
    """
    Calculates orthometric parameter "H3" in terms of mass and inclination.
    """

    return T_sun * m2 * stig_GR(sini)**3

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

