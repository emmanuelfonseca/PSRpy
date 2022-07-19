#/! /usr/bn/python

from astropy.coordinates import SkyCoord
from subprocess import call, Popen, PIPE
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import sys
import re

from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_name('serif')

parser = argparse.ArgumentParser(description='Generate a uniform, 2-D grid of TEMPO2 chi-squared values for different combinations of M2/COSI Shapiro parameters. Note for the future: work in an option to use M2/M1 (or M2/MTOT) combination instead.')

parser.add_argument('-f', nargs=1, action='store', dest='parfile', required=True, help='Input TEMPO parfile with M2/SINI parameters.')
parser.add_argument('-G', action='store_true', dest='useGLS', help='Use GLS fitting in tempo.')
parser.add_argument('-n', nargs=1, action='store', dest='Ngrid', default=[50], type=int, help='Number of grid points to generate. (Default: 50)')
parser.add_argument('-p', nargs=1, action='store', dest='pnum', default=[None], type=str, help='TEMPO pulse-number file. (Not required.)')
parser.add_argument('timfile', nargs=1, action='store', help='Input TEMPO2 .tim file containing TOAs.')
parser.add_argument('--RAJ', nargs=2, action='store', metavar=('min', 'max'), dest='RAJlimits', default=[-30, 30.], type=float, help='Lower/upper limits of position offsets along RAJ coordinate, in arcseconds. (Default: -30 to 30)')
parser.add_argument('--DECJ', nargs=2, action='store', metavar=('min', 'max'), dest='DECJlimits', default=[-30, 30.], type=float, help='Lower/upper limits of position offsets along DECJ coordinate, in arcseconds. (Default: -30 to 30)')
parser.add_argument('--tempo', action='store_true', dest='useTEMPO', default=False, help='If set, use TEMPO for gridding calculations.')
parser.add_argument('--qrfit', action='store_true', dest='useQRFIT', default=False, help='If set, use QR decomposition for TEMPO2 run.')

args = parser.parse_args()
intim = (args.timfile)[0]
inpar = (args.parfile)[0]
inpnum = (args.pnum)[0]
Ngrid = (args.Ngrid)[0]
ra_lo, ra_hi = args.RAJlimits
dec_lo, dec_hi = args.DECJlimits
useGLS = args.useGLS
useTEMPO = args.useTEMPO

ra_offsets = np.linspace(ra_lo, ra_hi, num=Ngrid)
dec_offsets = np.linspace(dec_lo, dec_hi, num=Ngrid)

# should not have to edit things from now on.
# but, you never know...

# set grid/tempo parameters.
chisq = np.zeros((Ngrid, Ngrid))
ra_bestfit = ''
dec_bestfit = ''

# get important parameters from TEMPO2 parfile.
for line in open(inpar, "r").read().splitlines():
    line = line.rstrip()
    if (re.search('PSR ', line)):
        line = line.split()
        obj = line[1]
    elif (re.search('PSRJ', line)):
        line = line.split()
        obj = line[1]
    elif (re.search('RAJ ', line)):
        line = line.split()
        ra_bestfit = line[1]
    elif (re.search('DECJ ', line)):
        line = line.split()
        dec_bestfit = line[1]

position_bestfit = SkyCoord(ra_bestfit, dec_bestfit, frame='icrs', unit=(u.hourangle, u.deg))

com1, com2 = [], []
chisq_bestfit = 0.
chisq_min = 1e6

if (useTEMPO):
    if useGLS:
        if (inpnum is not None):
            com1 = ['tempo', '-ni', inpnum, '-G', '-f', inpar, intim]
            com2 = ['tempo', '-ni', inpnum, '-G', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-ni', inpnum, '-G', '-f', './dummygrid2.par', intim]
        else:
            com1 = ['tempo', '-G', '-f', inpar, intim]
            com2 = ['tempo', '-G', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-G', '-f', './dummygrid2.par', intim]
    else:
        if (inpnum is not None):
            com1 = ['tempo', '-ni', inpnum, '-f', inpar, intim]
            com2 = ['tempo', '-ni', inpnum, '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-ni', inpnum, '-f', './dummygrid2.par', intim]
        else:
            com1 = ['tempo', '-f', inpar, intim]
            com2 = ['tempo', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-f', './dummygrid2.par', intim]
else:
    if useQRFIT:
        com1 = ['tempo2', '-qrfit', '-newpar', '-f', './dummygrid.par', intim]
        com2 = ['tempo2', '-qrfit', '-newpar', '-f', './dummygrid2.par', intim]
    else:
        com1 = ['tempo2', '-newpar', '-f', './dummygrid.par', intim]
        com2 = ['tempo2', '-newpar', '-f', './dummygrid2.par', intim]

print("Getting position grid...")

if useGLS:
    print("    * using GLS...")

current_ra_deg = 0.
current_dec_deg = 0.
c2 = 0
c3 = 0
chisq_min = 1e9
best_offset_ra = 0.
best_offset_dec = 0.

# loop over offsets; first in RAJ.
for current_offset_ra in ra_offsets:
    current_ra_deg = position_bestfit.ra.deg + (current_offset_ra / 3600.)
    c1 = 0

    # now loop over DECJ offsets.
    for current_offset_dec in dec_offsets:
        # put adjusted coordinates in appropriate string format.
        current_dec_deg = position_bestfit.dec.deg + (current_offset_dec / 3600.)
        current_position = SkyCoord(current_ra_deg, current_dec_deg, frame='icrs', unit=(u.deg, u.deg))
        current_ra_string = current_position.ra.to_string(unit=u.hourangle, precision=5, pad=True)
        current_ra_string = current_ra_string.replace("h", ":")
        current_ra_string = current_ra_string.replace("m", ":")
        current_ra_string = current_ra_string.replace("s", "")
        current_dec_string = current_position.dec.to_string(unit=u.degree, precision=5, pad=True)
        current_dec_string = current_dec_string.replace("d", ":")
        current_dec_string = current_dec_string.replace("m", ":")
        current_dec_string = current_dec_string.replace("s", "")

        fout = open('./dummygrid.par','w')

        # create new parfile with fixed M2/SINI (or M2/MTOT) values.
        for line in open(inpar, "r").read().splitlines():
            line = line.rstrip()

            if (re.search('RAJ ',line)):
                elem = line.split()
                fout.write("{0}               {1}  0\n".format('RAJ', current_ra_string))
            elif (re.search('DECJ ', line)):
                elem = line.split()
                fout.write("{0}               {1}  0\n".format('DECJ', current_dec_string))
            else:
                fout.write(line + "\n")

        fout.close()

        # now run tempo on new/dummy parfile, extract numbers.
        if (useTEMPO):

            p = Popen(com2, stdout=PIPE)
            out, err = p.communicate()
            try:
                m = re.search(r'Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', str(out))
                chisq[c1, c2] = float(m.group(3)) * float(m.group(2))
                if (chisq[c1, c2] < chisq_min):
                    chisq_min = chisq[c1, c2]
                    best_offset_ra = current_offset_ra
                    best_offset_dec = current_offset_dec
            except:
                chisq[c1, c2] = 1e6

        else:

            #call(com1, stdout=PIPE)
            p2 = Popen(com1, stdout=PIPE)
            out2, err2 = p2.communicate()
            p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
            out3, err3 = p3.communicate()
            line = out3.split()
            if (c1 == 0 and c2 == 0):
                nfree = float(line[2])
            chisq[c1,c2] = float(line[1]) * float(line[2])
            if (chisq[c1, c2] < chisq_min):
                chisq_min = chisq[c1, c2]
        c1 += 1
        c3 += 1
        perc = 100. * c3 / Ngrid**2
        sys.stdout.write('\rPercent completed: {0:.3f}% (min. chisq = {1:.3f}, corresponding [X, Y] = [{2:.3f}, {3:.3f}])'.format(perc, chisq_min, best_offset_ra, best_offset_dec))
        sys.stdout.flush()

    c2 += 1

print()

# create x/y grids for contour plots.
deltachi2 = chisq - np.min(chisq)
pdf2D = 0.5 * np.exp(-0.5 * deltachi2)

plt.pcolormesh(ra_offsets, dec_offsets, np.transpose(pdf2D), vmin=0, vmax=np.max(pdf2D), cmap="Blues")
plt.colorbar()
plt.xlabel(r'Offset in R. A. (arcsecond)', fontproperties=font, fontsize=15)
plt.ylabel(r'Offset in Dec. (arcescond)', fontproperties=font, fontsize=15)
plt.savefig("grid_ra_dec_{0}.png".format(obj), format='png')
plt.show()

GridDict = {}
GridDict['PSR'] = obj
GridDict['RAJ'] = ra_offsets
GridDict['DECJ']   = dec_offsets
GridDict['chisq'] = chisq
GridDict['RAJ_bestfit'] = position_bestfit.ra.deg
GridDict['DECJ_bestfit'] = position_bestfit.dec.deg

# save variables for future use.
pout = open("chisq_grid_position.pkl", 'wb')
pickle.dump(GridDict, pout)
pout.close()

# clean up shop.
q = Popen(['rm','dummygrid.par'],stdout=PIPE)
