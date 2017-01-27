"""
Listed below are the approximate orbital elements (and their first derivatives) 
of the major solar-system bodies. This list includes all 8 planets and Pluto. The 
values are taken from the document developed by E. M. Standish, and are valid 
for the time interval 1800 A.D. - 2050 A.D. Unless otherwise noted, all angular 
terms are given in units of degrees, and time derivatives are given with time units 
of centuries.

The following parameters are stored as nested Python dictionary entries with the 
following names for keys:

* sma     = semi-major axis [a.u.]
* sma-dot = time-derivative of sma [a.u. / cty]
* ecc     = eccentricity [  ]
* ecc-dot = time-derivative of ecc [    / cty]
* inc     = orbital inclination [deg]
* inc-dot = time-derivative of inc [deg / cty]
* mtl     = mean longitude [deg]
* mtl-dot = time-derivative of mtl [deg / cty]
* per     = perihelion longitude [deg]
* per-dot = time-derivative of per [deg / cty]
* asc     = longitude of ascending node [deg]
* asc-dot = time-derivative of asc
"""

ssb_planet_positions = {}
ssb_planet_positions['Earth'] = {
    "sma"     : 1.00000261,
    "sma-dot" : 0.00000562,
    "ecc"     : 0.01671123,
    "ecc-dot" : -0.00004392,
    "inc"     : -0.00001531, 
    "inc-dot" : -0.01294668,
    "mlt"     : 100.46457166,
    "mlt-dot" : 35999.37244981, 
    "per"     : 102.93768193, 
    "per-dot" : 0.32327364,
    "asc"     : 0.,
    "asc-dot" : 0.
}

