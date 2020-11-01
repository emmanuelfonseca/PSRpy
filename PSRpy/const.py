#! /usr/bin/python

import astropy.constants as ac
import astropy.units as u
import numpy as np

# physical quantities.
c            = ac.c                    # speed of light, m/s
G            = ac.G.value              # x Newton's constant, m^3 kg^-1 s^-2
G_err        = ac.G.uncertainty        # ... uncertainty
sigma_sb     = ac.sigma_sb             # Stefan-Boltzmann constant, W m^-2 K^-4
sigma_sb_err = ac.sigma_sb.uncertainty # ... uncertainty
k_B          = ac.k_B.value            # Boltzmann constant.
M_sun        = ac.M_sun.value          # mass of Sun, kg
M_sun_err    = ac.M_sun.uncertainty    # ... uncertainty
R_sun        = ac.R_sun.value          # radius of Sun, m
R_sun_err    = ac.R_sun.uncertainty    # ... uncertainty.
L_sun        = ac.L_sun.value          # luminosity of Sun, watts
T_sun        = G * M_sun / c**3        # Shapiro r for Sun, in units of seconds.
T_sun_err    = np.sqrt((G*M_sun_err/c**3)**2+(G_err*M_sun/c**3)**2)
au           = ac.au.value             # astronomical unit, m
pc           = ac.pc.value             # parsec, m
R0           = 8.34 * u.kpc            # distance to Galactic center, kpc.
R0_err       = 0.16 * u.kpc            # distance error, kpc
v0           = 240. * u.km / u.s      # circular speed of SSB, km/s
v0_err       = 8. * u.km / u.s        # circular speed error, km/s

# constants specific to the pulsar community.
k_DM = (1 / 0.000241) * u.s * u.MHz * u.MHz * u.cm * u.cm * u.cm / u.pc

# unit conversions.
d2r = np.pi / 180
r2d = 1 / d2r
spy = 86400 * 365.25
