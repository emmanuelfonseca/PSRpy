#! /usr/bin/python

from parfile import readpar, derivepar
from pkcorr import doppler
import numpy as np
import ddgr

par = readpar('./1909-3744.par')
der = derivepar(par)

b, l, mu = -0.3420093136832233, -0.004699546498681362, der.mu

m2, sini = par.M2, par.SINI
m1 = np.sqrt((m2*sini)**3/der.massfunc)-m2

# print stuff.
dop = par.PB*86400.*doppler(der.dist,b,l,mu)
print 'Contributions to PBDOT:'
print '    * Gal-accel:  {0:.5f}'.format(dop[0]*1e12)
print '    * Gal-relat:  {0:.5f}'.format(dop[1]*1e12)
print '    * Shklovksii: {0:.5f}'.format(dop[2]*1e12)
print 'PBDOT (GR):       {0:.5f}'.format(ddgr.pbdot(m1,m2,par.PB,par.E)*1e12)
print 'PBDOT (total):    {0:.5f}'.format((ddgr.pbdot(m1,m2,par.PB,par.E)+sum(dop))*1e12)
print 'PBDOT (observed): {0:.5f} +/- {1:.5f}'.format(par.PBDOT,par.PBDOTerr)
#print (ddgr.pbdot(m1,m2,par.PB,par.E)+sum(dop))*1e12

