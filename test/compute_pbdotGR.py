#! /usr/bin/python

from PSRpy.parfile import readpar, derivepar
from PSRpy.ddgr import doppler, PBDOT_GR
import numpy as np
import sys

path = '/Users/emmanuelfonseca/Documents/research/nanograv/binaryMSPs/parfiles/data_release/'

# compute for J0613-0200
parfile = path+'J0613-0200_NANOGrav_9yv1.gls.par'
par = readpar(parfile)
der = derivepar(par)

print "For J0613-0200: {0:.5f} (1e-12)".format(PBDOT_GR(1.57658, 0.15749, par.PB, der.E)*1e12)
