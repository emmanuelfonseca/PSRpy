#! /usr/bin/python

from PSRpy.parfile import readpar, derivepar, dmxpar, printpar
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
import sys

def pardict(parobj):
    outdict = {}
    for parameter in parobj.parorder:
        if (hasattr(parobj,parameter+'err')):
            outdict[parameter] = []
    return outdict

par = readpar('./J1949+3106.par')
der = derivepar(par)

print "(h3, stig) = ({0:.5f}, {1:.5f})".format(par.H3*Decimal(1e6),par.STIG)
print "(m2, sini) = ({0:.5f}, {1:.5f})".format(der.M2, der.SINI)
