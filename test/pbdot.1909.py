#! /usr/bin/python

from PSRpy.parfile import readpar, derivepar
from PSRpy.ddgr import doppler, PBDOT_GR
import numpy as np
import sys

par = readpar('./J1909.par')
der = derivepar(par)
b, l = -0.3420093136832233, -0.004699546498681362
mu, mf = float(der.mu), float(der.massfunc)
pb, PBD, PBDerr = float(par.PB), float(par.PBDOT), float(par.PBDOTerr)
ecc = float(der.E)
dist = float(der.dist)

muerr, disterr = float(der.muerr), float(der.disterr)

m2, sini = float(par.M2), float(par.SINI)
m1 = np.sqrt((m2 *sini)**3 / mf ) - m2

# print stuff.
dop, err = doppler(dist, disterr, b, l, mu, muerr)
dop, err = pb * 86400 * np.array([dop, err])
print 'Contributions to PBDOT (in units of 1e-12):'
print '    * Gal-potenial:   {0:.5f}'.format(dop[0]*1e12)
print '    * Gal-relative:   {0:.5f}'.format(dop[1]*1e12)
print '    * Shklovksii:     {0:.5f}'.format(dop[2]*1e12)
print 'PBDOT (GR, expected): {0:.5f}'.format(PBDOT_GR(m1,m2,pb,ecc)*1e12)
print 'PBDOT (total):        {0:.5f}'.format((PBDOT_GR(m1,m2,pb,ecc)+sum(dop))*1e12)
print 'PBDOT (observed):     {0:.5f} +/- {1:.5f}'.format(PBD, PBDerr)
