#! /usr/bin/python

### CONST: a python module containing constants and 
###        their uncertainties, when applicable.
### Emmanuel Fonseca, 24 January 2014.
import numpy as np

# physical quantities.
c       = 2.99792458e8   # speed of light, m/s
G       = 6.67428e-11    # Newton's constant, m^3 kg^-1 s^-2
Gerr    = 0.00067e-11    # ... uncertainty
sbc     = 5.670373e-8    # Stefan-Boltzmann constant, W m^-2 K^-4
sbcerr  = 0.000021e-8    # ... uncertainty
kb      = 1.38e-23
Msun    = 1.9884e30      # mass of Sun, kg
Msunerr = 0.0002e30      # ... uncertainty
Rsun    = 6.955e8        # radius of Sun, m
Teffsun = 5780           # temperature of Sun, kelvin
Lsun    = 3.856e26       # luminosity of Sun, watts
au      = 1.495978707e11 # astronomical unit, m
pc      = 3.08567758e16  # parsec, m
R0      = 8.34           # distance to Galactic center, kpc.
R0err   = 0.16           # distance error, kpc
v0      = 240.           # circular speed of SSB, km/s
v0err   = 8.             # circular speed error, km/s


# derived constants.
Tsun    = G*Msun/c**3    # 10^-6 s
Tsunerr = np.sqrt((G*Msunerr/c**3)**2+(Gerr*Msun/c**3)**2)
d2r     = np.pi/180.
spy     = 86400.*365.25
