#! /usr/bin/python

import matplotlib.pyplot as plt
from PSRpy.parfile import readpar
from re import search

inp = '/Users/emmanuelfonseca/Documents/research/1643-12/parfile/1643-1224.final.par'
par = readpar(inp)

dmx    = []
dmxerr = []
mjd    = []
mjdr1  = []
mjdr2  = []
f1     = []
f2     = []

for parameter in dir(par):
    if(search('DMX_',parameter) and not search('err',parameter)):
        line = parameter.split('_')
        dmx.append(getattr(par,'DMX_'+line[1]))
        dmxerr.append(getattr(par,'DMX_'+line[1]+'err'))
        mjd.append(getattr(par,'DMXEP_'+line[1]))
        mjdr1.append(getattr(par,'DMXR1_'+line[1]))
        mjdr2.append(getattr(par,'DMXR2_'+line[1]))
        f1.append(getattr(par,'DMXF1_'+line[1]))
        f2.append(getattr(par,'DMXF2_'+line[1]))

plt.plot(mjd,dmx,'r+')
plt.show()
