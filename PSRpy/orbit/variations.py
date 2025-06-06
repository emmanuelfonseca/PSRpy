#! /usr/bin/python

from ..const import c, G, M_sun, T_sun
import numpy as np

pi = np.pi

# compute the nominal GR parameters.

def gamma_GR(m1, m2, pb, e):
    """
    Calculates time-averaged gravitational redshift / time dilation parameter, 
    as expected from GR.
    """
    
    pb_in = pb * 86400
    A = e * (pb_in / 2 / pi)**(1./3.) * T_sun.value**(2./3.)
    return A * m2 * (m1 + 2 * m2) / (m1 + m2)**(4./3.)

def omdot_GR(m1, m2, pb, e, use_PK2=False, verbose=False):
    """
    Calculates periastron advance, as expected from GR. The default is to compute 
    the first-order PK term, but the second-order term can be computed if desired.
    """

    # first compute first-order term.
    pb_in = pb * 86400
    mtot = m1 + m2
    omdot = 3 * (pb_in / 2 / pi)**(-5./3.) * (T_sun.value * mtot)**(2./3.) / (1 - e**2)
    omdot_2PK = 0.

    # now compute and incorporate second-order term if desired.
    if use_PK2:
        x1 = (m1 / mtot)**2
        x2 = (m2 / mtot)**2
        x3 = m1 * m2 / mtot**2
        f0 = (39./4. * x1 + 27./4. * x2 + 15. * x3) / (1 - e**2) - \
             (13./4. * x1 + 1./4. * x2 + 13./3. * x3)
        b0 = (2 * pi * T_sun.value * mtot / pb_in) ** (1./3.)
        

        omdot_2PK = omdot * f0 * b0**2

    A = omdot * 180 / pi * 86400 * 365.25
    B = omdot_2PK * 180 / pi * 86400 * 365.25

    if verbose:
        print(f"OMDOT, first order: {A}")
        print(f"OMDOT, second order: {B}")
 
    return A + B

def pbdot_GR(m1, m2, pb, e, use_PK2=False, verbose=False):
    """
    Calculates orbital decay of binary period, as expected from GR.
    """

    pb_in = float(pb)*86400.
    nu = m1 * m2 / (m1 + m2) ** 2
    beta0 = (T_sun.value * (m1 + m2) * 2 * pi / pb_in) ** (1./3.)

    fe = (1 + (73. / 24. * e**2) + (37. / 96. * e**4)) / (1 - e**2) ** (7 / 2)
    ge = 0
    
    if use_PK2:
        ge = beta0 ** 2 / 336 / (1 - e ** 2) ** (9 / 2) * (
                 1273 + 16495 / 2 * e ** 2 + 
                 42231 / 8 * e ** 4 + 
                 3947 / 16 * e ** 6 - 
                 (
                     924 + 3381 * e ** 2 + 1659 / 4 * e ** 4 - 259 / 4 * e ** 6
                 ) * nu +
                 (
                     3297 * e ** 2 + 4221 * e ** 4 + 2331 / 8 * e ** 6
                 ) * (m1 - m2) / (m1 + m2)
             ) 

    #A = -192 * pi / 5 * (pb_in / 2 / pi)**(-5./3.) * fe * T_sun.value**(5./3.)
    A = -192 * pi / 5 * nu * (beta0 ** 5) * fe
    B = -192 * pi / 5 * nu * (beta0 ** 5) * ge

    if verbose:
        print(f"PBDOT, first order: {A}")
        print(f"PBDOT, second order: {B}")

    return A + B

def r_GR(m2):
    """
    Calculates Shapiro 'range' parameter, as expected from GR.
    """

    return T_sun * m2

def s_GR(m1, m2, pb, x):
    """
    Calculates Shapiro 'shape' parameter, as expected from GR / mass function. 
    Note that this is equal to the sine of the inclination angle.
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
    A = -64. / 5. * (2 * pi * T_sun.value / pb_in)**2 * fe
    return A * m1 * m2**2 / (m1 + m2)

def edot_GR(m1, m2, pb, e):
    """
    Calculates rate of circularization from GR.
    """

    nb = 2 * np.pi / pb / 86400
    eterm = (1 + (121./304.) * e**2) / (1 - e**2)**(5./2.)
    A = -304. / 15. * T_sun.value**(5./3.) * nb**(8./3.) * eterm
    return A * m1 * m2 / (m1 + m2)**(1./3.) * e

def dtheta_GR(m1, m2, pb, e):
    """
    Calculate the dtheta orbital-shape correction due to GR.
    """
    pb_in = pb * 86400
    nb = 2 * np.pi / pb_in
    A = (T_sun.value * nb)**(2./3.)
    return A * (7. * m1**2 / 2. + 6 * m1 * m2 + 2 * m2**2) / (m1 + m2)**(4./3.)

def precession_GR(m1, m2, pb, e):
    """
    Calculate rate of geodetic precession of the spin axis.
    """
    pb_in = pb * 86400
    A = (2 * np.pi / pb_in)**(5./3.) * T_sun.value**(2./3.) / (1 - e**2)
    return A * m2 * (4 * m1 + 3 * m2) / 2 / (m1 + m2)**(4./3.)


# compute effects and terms associated with Lense-Thirring precession.

def didt_LT(pb, e, m1, m2, moi, ps, pol_lambda, pol_eta):

    sl = np.sin(pol_lambda * pi / 180)
    ce = np.cos(pol_eta * pi / 180)
    mtot = (m1 + m2) * M_sun
    pb_in = pb * 86400
    moi_in = moi * 1e38
    sigma_p = G * (2 + 3 * m2 / m1 / 2) / c.value ** 2
    mean_motion = 2 * pi / pb_in
    sma_rel = (G * mtot / mean_motion ** 2) ** (1./3.)
    didt = sigma_p * moi_in * (2 * pi / ps) / (1 - e ** 2) ** (3./2.) / sma_rel ** 3 * sl * ce

    return didt

def omdot_LT(pb, e, m1, m2, moi, ps, pol_lambda, pol_eta, incl):
    
    pb_in = pb * 86400
    mtot = m1 + m2
    x1 = m1 / mtot
    x2 = m2 / mtot
    moi_in = moi * 1e38
    beta_spin = G * moi_in * (2 * pi / ps) / (m1 * T_sun.value) ** 2 / c.value ** 5
    mean_motion = 2 * pi / pb_in
    beta_orbit = (mean_motion * mtot * T_sun.value) ** (1./3.)

    # now compute the term that depends on geometry.
    sine = np.sin(pol_eta * pi / 180)
    cose = np.cos(pol_eta * pi / 180)
    sini = np.sin(incl * pi / 180)
    cosi = np.cos(incl * pi / 180)
    sinl = np.sin(pol_lambda * pi / 180)
    cosl = np.cos(pol_lambda * pi / 180)

    cos_delta = -sini * sinl * sine + cosl * cosi
    g_spin = x1 * (4 * x1 + 3 * x2) / 6 / np.sqrt(1 - e ** 2) / sini ** 2
    g_spin *= ((3 * sini ** 2 - 1) * cos_delta + cosi * cosl)

    # finally, compute OMDOT due to LT precession.
    omdot = omdot_GR(m1, m2, pb, e)
    omdot_LT = omdot * beta_orbit * beta_spin * g_spin 

    return omdot_LT 

def xdot_LT(pb, e, m1, m2, moi, ps, pol_lambda, pol_eta, a1, sini):

    cosi = np.sqrt(1 - sini ** 2)
    coti = cosi / sini
    didt = didt_LT(pb, e, m1, m2, moi, ps, pol_lambda, pol_eta)
    xdot = a1 * coti * didt

    return xdot


# compute terms from loss of spin/mass energy.
def pbdot_mdot(pb, m1, m2, moi, ps, psdot):

    pb_in = pb * 86400
    mtot = (m1 + m2) * M_sun
    moi_in = moi * 1e38
    pbdot = 8 * pi ** 2 * moi_in * psdot * pb_in / mtot / ps ** 3 / c.value ** 2
  
    return pbdot  

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

    return T_sun.value * m2 * stig_GR(sini)**3

# compute terms for geodetic precession.
def precession_geod(pb, e, m1, m2):
    
    pb_in = pb * 86400
    mtot = m1 + m2
    mean_motion = 2 * pi / pb_in
    A = T_sun.value ** (2. / 3) * mean_motion ** (5. / 3) / 2 / (1 - e ** 2)
    prec = A * m2 * (4 * m1 + 3 * m2) / mtot ** (4. / 3)
    prec *= (180 / pi * 86400 * 365.25)
    
    return prec

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

def xdot_aberration(pb, e, m1, m2, ps, pol_lambda, pol_eta, a1, sini):

    prec_in = precession_geod(pb, e, m1, m2) * (pi / 180) / 365.25 / 86400
    cosi  = np.sqrt(1 - sini ** 2)
    sine  = np.sin(pol_eta * pi / 180)
    cose  = np.cos(pol_eta * pi / 180)
    sin2e = np.sin(2 * pol_eta * pi / 180)
    cos2e = np.cos(2 * pol_eta * pi / 180)
    sinl  = np.sin(pol_lambda * pi / 180)
    cosl  = np.cos(pol_lambda * pi / 180)
    

    d_eps_dt = -ps / (pb * 86400) / sini / np.sqrt(1 - e ** 2)
    d_eps_dt *= prec_in / sinl ** 2
    d_eps_dt *= (sini * cosl * sin2e + cosi * sinl * cose)
    xdot = a1 * d_eps_dt

    return xdot

# compute terms due to change in orientation from proper motion.

def xdot_pm(pmra, pmdec, longitude_node_asc, a1, incl):
    pm = np.sqrt(pmra ** 2 + pmdec ** 2)
    pm *= (pi / 180 / 1000 / 3600 / 365.25 / 86400)
    sini = np.sin(incl * pi / 180)
    cosi = np.cos(incl * pi / 180)
    coti = cosi / sini
    posang_pm = np.arctan2(pmdec, pmra) * 180 / pi 
    diff = (posang_pm - longitude_node_asc) * pi / 180
    xdot = a1 * pm * coti * np.sin(diff)    

    return xdot


