#! /usr/bin/python

from subprocess import call, Popen, PIPE
from PSRpy.const import T_sun
from os.path import isfile
import PSRpy.orbit.variations as orbvar
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import sys
import re

from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_name('serif')

parser = argparse.ArgumentParser(description=
    "Generates a uniform grid of TEMPO/TEMPO2 chi-squared values for different " +  
    "combinations of M2/COSI parameters of the Shapiro time delay. If certain " + 
    "command-line options are chosen, a three-dimensional grid will be computed, " + 
    "with the third axis being either the pulsar distance or ascending-node longitude."
)

parser.add_argument(
    "timfile", action="store", 
    help="ASCII file containing TOAs that are compatible with TEMPO or TEMPO2."
)

parser.add_argument(
    "-f", action='store', dest="parfile", required=True, 
    help="Input TEMPO or TEMPO2 parfile with M2/SINI parameters."
)

parser.add_argument(
    "-G", action="store_true", dest="useGLS", help="Use GLS fitting in TEMPO."
)

parser.add_argument(
    "-n", action="store", dest="Ngrid", default=50, type=int, 
    help="Number of grid points to generate. (Default: 50)"
)

parser.add_argument(
    "-p", action="store", dest="pnum", default=None, type=str, help="TEMPO pulse-number file."
)

parser.add_argument(
    "--M2", nargs=2, action="store", metavar=("min", "max"), dest="M2limits", default=[0.01, 1.], 
    type=float, help="Limits of grid along M2 coordinate. (Default: 0.01 to 1)"
)

parser.add_argument(
    "--COSI", nargs=2, action="store", metavar=("min", "max"), dest="COSIlimits", 
    default=[0., 0.99], type=float, 
    help="Limits of grid along COSI coordinate. (Default: 0 to 0.99)"
)

parser.add_argument(
    "--H3", nargs=2, action="store", metavar=("min", "max"), dest="H3limits", 
    default=[0.01, 1.0], type=float, 
    help="Limits of grid along H3 coordinate. (Default: 0.01 to 1.0)"
)

parser.add_argument(
    "--H4", nargs=2, action="store", metavar=("min", "max"), dest="H4limits", 
    default=[0.01, 1.0], type=float, 
    help="Limits of grid along H4 coordinate. (Default: 0.01 to 1.0)"
)

parser.add_argument(
    "--STIG", nargs=2, action="store", metavar=("min", "max"), dest="STIGlimits", 
    default=[0.01, 0.99], type=float, 
    help="Limits of grid along STIG coordinate. (Default: 0.01 to 0.99)"
)

parser.add_argument(
    "--THETA", nargs=2, action="store", metavar=("min", "max"), dest="THETAlimits", 
    default=[0., 0.], type=float, 
    help="Limits of grid along THETA (kinematic XDOT) coordinate. (Default: 0 to 360 degrees.)"
)

parser.add_argument("--M1", nargs=2, action="store", metavar=("min", "max"), dest="M1limits", 
    default=[0.1, 10.], type=float, 
    help="Limits of grid along M1 coordinate. (Default: 0 to 10 solar masses.)"
)

parser.add_argument(
    "--MTOT", nargs=2, action="store", metavar=("min", "max"), dest="MTOTlimits", 
    default=[0.1, 10.], type=float, 
    help="Limits of grid of values for MTOT. (Default: 0 to 10 solar masses.)"
)

parser.add_argument(
    "--PX", nargs=2, action="store", metavar=("min", "max"), dest="PXlimits", 
    default=[0., 0.], type=float, 
    help="Limits of grid along PX coordinate. (Default: 0.01 to 10 mas.)"
)

parser.add_argument(
    "--XOMDOT", nargs=2, action="store", metavar=("min", "max"), dest="XOMDOTlimits", 
    default=[0., 0.], type=float, 
    help="Lower/upper limits of grid of values for PX. (Default: 0.01 to 10 mas.)"
)

parser.add_argument(
    "--H3STIG", action="store_true", dest="gridH3STIG", 
    help="Grid over the H3/STIG parameters. (Requires a solution that uses these parameters.)"
)

parser.add_argument(
    "--M2MTOT", action="store_true", dest="gridM2MTOT", 
    help="Grid of the M2/MTOT parameters instead of M2/COSI."
)

parser.add_argument(
    "--OMDOT", action="store_true", dest="fixOMDOT", 
    help="Compute and fix GR component of OMDOT for each M2/COSI coordinate. " + 
         "(This option only works if OMDOT is set in the input parfile.)"
)

parser.add_argument("--PBDOT", action="store_true", dest="fixPBDOT", 
    help="Compute and fix GR component of PBDOT for each M2/COSI coordinate. " + 
         "(This option only works if PBDOT is set in the input parfile.)'"
)

parser.add_argument(
    "--GAMMA", action="store_true", dest="fixGAMMA", 
    help="Compute and fix GR component of GAMMA for each M2/COSI coordinate. " + 
         "(This option only works if GAMMA is set in the input parfile.)"
)

parser.add_argument(
    "--DDGR", action="store_true", dest="gridDDGR", 
    help="Grid over M2/MTOT using the DDGR binary model; " + 
         "currently floats all other DDGR parameters (e.g. XOMDOT, XPBDOT)."
)

parser.add_argument(
    "--H3H4", action="store_true", dest="gridH3H4", 
    help="Grid over H3/H4 using the ELL1H binary model."
)

parser.add_argument(
    "--DDS", action="store_true", dest="gridDDS", 
    help="Grid over M2/COSI using the DDS binary model."
)

parser.add_argument(
    "--DDK", action="store_true", dest="gridDDK", 
    help="Grid over M2/COSI/KOM using the DDK binary model."
)

parser.add_argument(
    "--tempo", action="store_true", dest="useTEMPO", 
    help="Perform grid using TEMPO instead."
)

parser.add_argument(
    "--qrfit", action="store_true", dest="useQRFIT", 
    help="Use QR matrix decomposition in TEMPO2 instead of Cholesky method."
)

parser.add_argument(
    "--M1M2", action="store_true", dest="grid_m1m2", help='Uniformly grid over M1/M2.'
)

parser.add_argument("--PK2", action="store_true", dest="use_PK2", 
    help="Compute second order post-Keplerian version of OMODT."
)

parser.add_argument(
    "--EXP1", action="store_true", dest="exp_OMDOT", 
    help="Skip TEMPO run if OMDOT is outside of 30-sigma range. Experimental feature."
)

parser.add_argument(
    "--fitM2", action="store_true", dest="fitM2", help="Force M2 to be fit during grid."
)

parser.add_argument(
    "--fitSINI", action="store_true", dest="fitSINI", help="Force SINI to be fit during grid."
)

# collect and retrieve command-line arguments.
args = parser.parse_args()
intim = args.timfile
inpar = args.parfile
inpnum = args.pnum
Ngrid = args.Ngrid
m2_lo, m2_hi = args.M2limits
cosi_lo, cosi_hi = args.COSIlimits
h3_lo, h3_hi = args.H3limits
h4_lo, h4_hi = args.H4limits
stig_lo, stig_hi = args.STIGlimits
theta_lo, theta_hi = args.THETAlimits
m1_lo, m1_hi = args.M1limits
mtot_lo, mtot_hi = args.MTOTlimits
px_lo, px_hi = args.PXlimits
xomdot_lo, xomdot_hi = args.XOMDOTlimits
fixOMDOT = args.fixOMDOT
fixPBDOT = args.fixPBDOT
fixGAMMA = args.fixGAMMA
fixXOMDOT = False
fixXDOT = False
fixPX = False
gridDDGR = args.gridDDGR
gridDDS = args.gridDDS
gridDDK = args.gridDDK
gridH3H4 = args.gridH3H4
gridH3STIG = args.gridH3STIG
gridM1M2 = args.grid_m1m2
gridM2MTOT = args.gridM2MTOT
useTEMPO = args.useTEMPO
useQRFIT = args.useQRFIT
useGLS = args.useGLS
use_PK2 = args.use_PK2
expOMDOT = args.exp_OMDOT
fitM2 = args.fitM2
fitSINI = args.fitSINI
grid_m2cosi = True

# decide which type of grid to generate.
if (gridDDGR or gridM2MTOT):
    grid_m2cosi = False

if (any([px_lo, px_hi]) != 0.):
    fixPX = True

if (any([theta_lo, theta_hi]) != 0.):
    fixXDOT = True

if (any([xomdot_lo, xomdot_hi]) != 0.):
    fixXOMDOT = True

if (gridM1M2 or gridH3STIG or gridH3H4):
    grid_m2cosi = False

# define possible grid params here.
# TODO: make this step more efficient by only initializing 
# X, Y, and (if applicable) Z arrays, and overloading those 
# values with arrays for the appropriate parameters.
cosi = np.linspace(cosi_lo, cosi_hi, num=Ngrid)
h3 = np.linspace(h3_lo, h3_hi, num=Ngrid)
h4 = np.linspace(h4_lo, h4_hi, num=Ngrid)
m1 = np.linspace(m1_lo, m1_hi, num=Ngrid)
m2 = np.linspace(m2_lo, m2_hi, num=Ngrid)
mtot = np.linspace(mtot_lo, mtot_hi, num=Ngrid)
px = np.linspace(px_lo, px_hi, num=Ngrid)
stig = np.linspace(stig_lo, stig_hi, num=Ngrid)
theta = np.linspace(theta_lo, theta_hi, num=Ngrid)
xomdot = np.linspace(xomdot_lo, xomdot_hi, num=Ngrid)

# set grid/tempo parameters.
tolerance = 30.
chisq = np.zeros((Ngrid, Ngrid))
chisq_3D = np.zeros((Ngrid, Ngrid, Ngrid))
sini_bestfit, m2_bestfit, mtot_bestfit = 0, 0, 0
h3_bestfit, h4_bestfit, stig_bestfit = 0., 0., 0.
Pb, A1, E = 0, 0, 0
xdot, pbdot, omdot, gamma = 0., 0., 0., 0.
pmlambda, pmbeta = 0., 0.
pmra, pmdec = 0., 0.
pm, px_bestfit = 0., 0.
gal_l, gal_b = 0., 0.
nfree, nharm = 0, 0
omdoterr = 0.
obj = ""
binary_model = ""

# get important parameters from TEMPO2 parfile.
for line in open(inpar, "r").read().splitlines():
    line = line.rstrip()
    if (re.search('PSR ', line)):
        line = line.split()
        obj = line[1]
    elif (re.search('PSRJ', line)):
        line = line.split()
        obj = line[1]
    elif (re.search('PMRA ', line)):
        line = line.split()
        pmra = np.float(line[1])
    elif (re.search('PMDEC ', line)):
        line = line.split()
        pmdec = np.float(line[1])
    elif (re.search('PMLAMBDA ', line)):
        line = line.split()
        pmlambda = np.float(line[1])
    elif (re.search('PMBETA ', line)):
        line = line.split()
        pmbeta = np.float(line[1])
    elif (re.search('H3 ', line)):
        line = line.split()
        h3_bestfit   = np.float(line[1])
    elif (re.search('H4 ', line)):
        line = line.split()
        h4_bestfit   = np.float(line[1])
    elif (re.search('M2 ', line)):
        line = line.split()
        m2_bestfit   = np.float(line[1])
    elif (re.search('STIG ', line)):
        line = line.split()
        stig_bestfit = np.float(line[1])
    elif (re.search('SINI ', line)):
        line = line.split()
        sini_bestfit = np.float(line[1])
    elif (re.search('MTOT ', line)):
        line = line.split()
        mtot_bestfit = np.float(line[1])
    elif (re.search('PB ', line)):
        line = line.split()
        Pb = np.float(line[1])
    elif (re.search('A1 ', line)):
        line = line.split()
        A1 = np.float(line[1])
    elif (re.search('ECC ', line)):
        line = line.split()
        E = np.float(line[1])
    elif (line[0] == 'E' and line[1] == ' '):
        line = line.split()
        E = np.float(line[1])
    elif (re.search('OMDOT ', line)):
        line = line.split()
        omdot = np.float(line[1])
        if (len(line) > 2):
            omdoterr = np.float(line[3])
    elif (re.search('GAMMA ', line)):
        line = line.split()
        gamma = np.float(line[1])
    elif (re.search('XDOT ', line)):
        line = line.split()
        xdot = np.float(line[1])
    elif (re.search('BINARY ', line)):
        binary_model = line.split()[1]

massfunc = A1**3 * (2 * np.pi / Pb / 86400.)**2 / T_sun.value

if (pmra !=0 and pmdec != 0):
    pm = np.sqrt(pmra**2 + pmdec**2)
else:
    pm = np.sqrt(pmlambda**2 + pmbeta**2)

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

print("Getting Shapiro-delay grid...")

if useGLS:
    print("    * using GLS...")

if fixOMDOT:
    if (omdot != 0.):
        print("    * fixing OMDOT to GR value (best-fit OMDOT = {0:.5f} deg/yr)".format(omdot))
    else:
        print("    * !!! No OMDOT in original parfile; proceed with caution !!!")
    #   sys.exit("Ooops! It does not look like OMDOT is set in this parfile.")

if fixPBDOT:
    if (hasattr(par, 'PBDOT') and fixPX):
        print("    * fixing PBDOT to GR+kinematic value (best-fit PBDOT = {0:.5f} * 1e-12)".format(pbdot))
    elif hasattr(par, 'PBDOT'):
        print("    * fixing PBDOT to GR value (best-fit PBDOT = {0:.5f} * 1e-12)".format(pbdot))
    else:
        sys.exit("Ooops! It does not look like PBDOT is set in this parfile.")

if fixGAMMA:
    if (gamma != 0.):
        print("    * fixing GAMMA to GR value (best-fit GAMMA = {0:.5f} ms)".format(gamma * 1000.))
    else:
        sys.exit("Ooops! It does not look like GAMMA is set in this parfile.")

if fixXDOT:
    print("Fixing XDOT (best-fit = {0:.5f})...".format(xdot))
    cosi_bestfit = np.sqrt(1 - sini_bestfit**2)
    coti_bestfit = cosi_bestfit / sini_bestfit
    xdot_max = pm * A1 * coti_bestfit * 1e12
    sinTHETA = xdot / xdot_max
    print("    * predicted max. XDOT = {0:.5f}...".format(xdot_max))
    print("    * predicted THETA     = {0:.5f}...".format(np.arcsin(sinTHETA) * 180 / np.pi))

if gridH3H4:
    print("    * orthometric (approximate) timing solution...")
    print("        - best-fit H3: {0:.3f} microsec".format(h3_bestfit * 1e6))
    print("        - best-fit H4: {0:.3f} microsec".format(h4_bestfit * 1e6))

if gridH3STIG:
    print("    * orthometric (full) timing solution...")
    print("        - best-fit H3:   {0:.3f} microsec".format(h3_bestfit * 1e6))
    print("        - best-fit STIG: {0:.3f}".format(stig_bestfit))
    sini_bestfit = 2 * stig_bestfit / (1 + stig_bestfit**2)
    m2_bestfit = h3_bestfit / stig_bestfit**3 / T_sun
    print("    * derived M2:    {0:.3f} solar masses".format(m2_bestfit))
    print("    * derived SINI:  {0:.3f}".format(sini_bestfit))

if (gridDDGR):
    print("Fixing MTOT (best-fit = {0:.5f})...".format(mtot_bestfit))

c2, c3 = 0, 0
mtot_parfile = 0.
x, y = 0, 0

# figure out which grid to make, select appropriate axes.
if grid_m2cosi:
    x = cosi
    y = m2

elif gridH3H4:
    x = h4
    y = h3

elif gridH3STIG:
    x = stig
    y = h3

elif (gridDDGR or gridM2MTOT):
    x = mtot
    y = m2

elif gridM1M2:
    x = m1
    y = m2

else:
    sys.exit("It's unclear which grid you want to generate!")

gamma_new = 0.
pbdot_new = 0.
xdot_new = 0.
omdot_new = 0.
omdot_full = 0.
sini_elem = 0.
cosi_elem = 0.
shapmax_elem = 0.
m1_elem = 0.
m2_elem = 0.
mtot_elem = 0.
count_3D = 0

# now perform loops over the axes.
for x_elem in x:

    # declare appropriate variables, depending on desired grid.
    if grid_m2cosi:
        cosi_elem = x_elem
        sini_elem = np.sqrt(1 - cosi_elem**2)
        shapmax_elem = -np.log(1 - sini_elem)

    elif (gridDDGR or gridM2MTOT):
        mtot_elem = x_elem

    elif gridM1M2:
        m1_elem = x_elem

    c1 = 0

    # loop over y coordinate.
    for y_elem in y:

        if (grid_m2cosi or gridDDGR or gridM1M2 or gridM2MTOT):
            m2_elem = y_elem

            if grid_m2cosi:
                m1_elem = np.sqrt((m2_elem * sini_elem)**3 / massfunc) - m2_elem

            elif (gridDDGR or gridM1M2 or gridM2MTOT):
                if (gridDDGR or gridM2MTOT):
                    m1_elem = mtot_elem - m2_elem

                sini_elem = (massfunc * (m1_elem + m2_elem)**2)**(1./3.) / m2_elem

        if (fixXDOT or fixPX):
            pass
        else:
            perc = 100.*(c3+1)/Ngrid**2
            sys.stdout.write('\rPercent completed: {0:.3f}% (min. chisq = {1:.3f})'.format(perc, chisq_min))
            sys.stdout.flush()

        fout = open('./dummygrid.par','w')

        # create new parfile with fixed M2/SINI (or M2/MTOT) values.
        for line in open(inpar, "r").read().splitlines():
            line = line.rstrip()

            # if DDGR solution, grid over different, DDGR-specific parameters.
            if gridDDGR:

                # uniform grid in M2/COSI? Or M1/M2? 
                #if grid_m2cosi:
                #    mtot_parfile = np.sqrt((m2_elem * sini_elem)**3 / massfunc)
                #else:
                #    mtot_parfile = m1_elem + m2_elem

                if (re.search('M2 ', line)):
                    line = line.split()
                    fout.write("{0}               {1:.12f}  0\n".format(line[0], m2_elem))
                elif (re.search('MTOT ', line)):
                    line = line.split()
                    fout.write("{0}               {1:.12f}  0\n".format(line[0], mtot_elem))
                elif ('GAMMA ' in line):
                    pass
                elif ('SINI ' in line):
                    pass
                elif ('OMDOT ' in line):
                    pass
                elif ('DR ' in line):
                    pass
                elif ('DTHETA ' in line):
                    pass
                elif ('XPBDOT' not in line and 'PBDOT' in line):
                    pass
                else:
                    fout.write(line+"\n")

            elif (gridDDS):
                if (re.search('M2 ', line)):
                    line = line.split()
                    fout.write("{0}               {1:.12f}  0\n".format(line[0], m2_elem))
                elif (re.search('SHAPMAX ', line)):
                    line = line.split()
                    fout.write("{0}               {1:.12f}  0\n".format(line[0], shapmax_elem))
                else:
                    fout.write(line+"\n")

            elif gridH3H4:
                if (re.search('H3 ', line)):
                    line = line.split()
                    fout.write("{0}       {1:.12f}  0\n".format(line[0], y_elem * 1e-6))
                elif (re.search('H4 ', line)):
                    line = line.split()
                    fout.write("{0}       {1:.12f}  0\n".format(line[0], x_elem * 1e-6))
                else:
                    fout.write(line+"\n")

            # else, if orthometric solution, grid over appropriate orthometric parameters.
            elif gridH3STIG:

                h3_elem = y_elem * 1e-6
                stig_elem = x_elem
                m2_elem = h3_elem / stig_elem**3 / T_sun
                sini_elem = 2 * stig_elem / (1 + stig_elem**2)
                m1_elem = np.sqrt((m2_elem * sini_elem)**3 / massfunc) - m2_elem
                elem = line.split()

                if (re.search('H3 ', line)):
                    line = line.split()
                    fout.write("{0}       {1:.12f}  0\n".format(line[0], y_elem * 1e-6))
                elif (re.search('VARSIGMA ', line)):
                    line = line.split()
                    fout.write("{0}       {1:.8f}  0\n".format(line[0], x_elem))
                elif (re.search('STIG ', line)):
                    line = line.split()
                    fout.write("{0}       {1:.8f}  0\n".format(line[0], x_elem))
                elif (re.search('OMDOT',line)):
                    if (fixOMDOT and fixXDOT):
                        m1_elem = np.sqrt((m2_elem * sini_elem)**3 / massfunc) - m2_elem
                        #omdot_new = ddgr.omdot_GR(m1_elem, m2_elem, Pb, E)
                    elif fixOMDOT:
                        pass
                        #omdot_new = ddgr.omdot_GR(m1_elem, m2_elem, Pb, E)
                        #fout.write("{0}               {1:.12f}  0\n".format(elem[0], omdot_new))
                    else:
                        fout.write(line+"\n")
                elif (re.search('PBDOT',line)):
                    if (fixPBDOT and fixPX):
                        #pbdot_new = ddgr.pbdot_GR(m1_elem, m2_elem, Pb, E)
                        pass
                    elif fixPBDOT:
                        m1_elem = np.sqrt((m2_elem * sini_elem)**3 / massfunc) - m2_elem
                        #pbdot_new = ddgr.pbdot_GR(m1_elem, m2_elem, Pb, E)
                        #fout.write("{0}               {1:.8f}  0\n".format(elem[0], pbdot_new * 1e12))
                    else:
                        fout.write(line+"\n")
                elif (re.search('GAMMA',line)):
                    if fixGAMMA:
                        #gamma_new = ddgr.gamma_GR(m1_elem, m2_elem, Pb, E)   
                        pass
                        #fout.write("{0}               {1:.12f}  0\n".format(elem[0], gamma_new))
                    else:
                        fout.write(line+"\n")
                elif (re.search('XDOT',line)):
                    if fixXDOT:
                        pass
                elif (re.search('PX',line)):
                    if fixPX:
                        pass
                    else:
                        fout.write(line+"\n")
                else:
                    fout.write(line+"\n")

            # otherwise, assumes DD/ELL1/BT/DDK parameters.
            else:

                if (re.search('SINI ',line)):
                    elem = line.split()
                    fout.write("{0}               {1:.8f}  0\n".format('SINI', sini_elem))
                elif (re.search('KOM ', line) and fixXDOT):
                    pass
                elif (re.search('KOM ', line) and gridDDK):
                    pass
                elif (re.search('KOM ', line)):
                    fout.write("{0}\n".format(line))
                elif (re.search('KIN ',line)):
                    incl_in = np.arccos(cosi_elem) * 180 / np.pi
                    fout.write("{0}               {1:.8f}  0\n".format('KIN',incl_in))
                elif (re.search('M2 ',line) and not re.search('DM2', line)):
                    if fitM2:
                        fout.write(line)
                    else:
                        elem = line.split()
                        fout.write("{0}                 {1:.8f}  0\n".format(elem[0], m2_elem))
                elif (re.search('OMDOT ',line)):
                    if ((fixOMDOT and fixXOMDOT) or (fixOMDOT and fixXDOT)):
                        omdot_new = orbvar.omdot_GR(m1_elem, m2_elem, Pb, E, use_PK2=use_PK2)
                    elif fixOMDOT:
                        omdot_new = orbvar.omdot_GR(m1_elem, m2_elem, Pb, E, use_PK2=use_PK2)
                        fout.write("{0}               {1:.12f}  0\n".format('OMDOT', omdot_new))
                    else:
                        fout.write(line+"\n")
                #elif (re.search('PBDOT',line)):
                #    if (fixPBDOT and fixPX):
                #        pbdot_new = ddgr.pbdot_GR(m1_elem, m2_elem, Pb, E)
                #    elif fixPBDOT:
                #        pbdot_new = ddgr.pbdot_GR(m1_elem, m2_elem, Pb, E)
                #        fout.write("{0}               {1:.8f}  0\n".format('PBDOT', pbdot_new * 1e12))
                #    else:
                #        fout.write(line+"\n")
                elif (re.search('GAMMA',line)):
                    if fixGAMMA:
                        gamma_new = orbvar.gamma_GR(m1_elem, m2_elem, Pb, E).value
                        fout.write("{0}               {1:.12f}  0\n".format('GAMMA', gamma_new))
                    else:
                        fout.write(line+"\n")
                elif (re.search('XDOT',line)):
                    if (fixXDOT):
                        # if 3-D grid is wanted, then omit the 'XDOT' line from parfile.
                        # this is because you need to set KIN/KOM params instead.
                        pass
                    else:
                        fout.write(line+"\n")
                elif (re.search('PX',line)):
                    if fixPX:
                        pass
                    else:
                        fout.write(line+"\n")
                else:
                    fout.write(line+"\n")

        # if 3D grid is wanted, then loop through dummy parfile injecting different values of KOM or PX.
        # otherwise, close filehandle and get chi-squared values.
        #if (omdot == 0. and fixOMDOT):
        #    omdot_new = ddgr.omdot_GR(m1_elem, m2_elem, Pb, E)
        #    fout.write("{0}             {1:.12f}  0\n".format("OMDOT", omdot_new))
        fout.close()

        if (gridDDK or fixXDOT):

            c4 = 0

            for theta_elem in theta:

                perc = 100. * (count_3D + 1)/Ngrid**3
                sys.stdout.write('\rPercent completed: {0:.3f}% (min. chisq = {1:.3f})'.format(perc, chisq_min))
                sys.stdout.flush()
                c3 += 1
                call(['cp', 'dummygrid.par', 'dummygrid2.par'])
                fout2 = open('./dummygrid2.par', 'a')

                if (fixXDOT and not gridDDK):
                    coti_elem = np.sqrt(1 - sini_elem**2) / sini_elem
                    xdot_full = A1 * pm * coti_elem * np.sin(theta_elem * np.pi / 180) / 1000. / 3600. 
                    xdot_full *= (np.pi / 180)
                    xdot_full /= (365.25 * 86400)
                    omdot_sec = pm / sini_elem * np.cos(theta_elem * np.pi / 180) / 1000. / 3600.
                    omdot_full = omdot_new + omdot_sec

                    # write out XDOT in correct units if using TEMPO or TEMPO2.
                    if useTEMPO:
                        fout2.write("{0}                 {1:.8f}   0\n".format('XDOT', xdot_full * 1e12))
                    else:
                        fout2.write("{0}                 {1:.8e}   0\n".format('XDOT', xdot_full))

                    # if model is not ELL1, write out OMDOT.
                    if (binary_model != "ELL1"):
                        fout2.write("{0}                {1:.8f}   0\n".format('OMDOT', omdot_full))

                # else, assume DDK model and set KOM.
                else:
                    fout2.write("{0}                 {1:.8f}   0\n".format('KOM', theta_elem))

                fout2.close()

                # now run tempo on new/dummy parfile, extract numbers.
                if useTEMPO:

                    p3 = Popen(com3, stdout=PIPE)
                    out3, err3 = p3.communicate()
                    try:
                        m = re.search(r'Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', out3.decode('utf-8'))
                        chisq_3D[c1, c2, c4] = np.float(m.group(3)) * np.float(m.group(2))
                        if (chisq_3D[c1, c2, c4] < chisq_min):
                            chisq_min = chisq_3D[c1, c2, c4]

                    except:
                        chisq_3D[c1,c2,c4] = 1e6

                    c4 += 1

                else:
                    call(com2,stdout=PIPE)
                    p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
                    out3, err3 = p3.communicate()
                    line = out3.split()
                    if (c1 == 0 and c2 == 0 and c4 == 0):
                        nfree = np.float(line[2])
                    chisq_3D[c1,c2, c4] = np.float(line[1]) * np.float(line[2])
                    c4 += 1

                count_3D += 1 

        elif (fixPX or fixXOMDOT):

            c4 = 0
            z = 0

            if fixPX:
                z = px
            elif fixXOMDOT:
                z = xomdot

            for z_elem in z:

                perc = 100. * (count_3D + 1) / Ngrid**3
                sys.stdout.write('\rPercent completed: {0:.3f}% (min. chisq = {1:.3f})'.format(perc, chisq_min))
                sys.stdout.flush()
                c3 += 1
                call(['cp', 'dummygrid.par', 'dummygrid2.par'])

                if fixPX:

                    fout2 = open('./dummygrid2.par', 'a')
                    fout2.write("{0}                 {1:.8f}   0\n".format('PX', z_elem))
                    
                    if (fixPBDOT):
                        #dop_comps, err = doppler(1 / z_elem, 0.2, gal_b, gal_l, pm, 0.1)
                        pbdot_full = pbdot_new #+ np.sum(dop_comps) * Pb * 86400
                        fout2.write("{0}                 {1:.8f}   0\n".format('PBDOT', pbdot_full * 1e12))

                    fout2.close()

                elif fixXOMDOT:

                    omdot_full = omdot_new + z_elem * 1e-6

                    fout2 = open('./dummygrid2.par', 'a')
                    fout2.write("{0}                 {1:.8f}   0\n".format('OMDOT', omdot_full))
                    fout2.close()


                # now run tempo on new/dummy parfile, extract numbers.

                if (fixPBDOT and expOMDOT and np.fabs(pbdot_full * 1e12 - pbdot) / pbdoterr > tolerance):
                    
                    chisq_3D[c1, c2, c4] = 1e6
                    count_3D += 1
                    c4 += 1

                elif (fixGAMMA and expOMDOT and np.fabs(gamma_new - gamma) / gammaerr > tolerance):

                    chisq_3D[c1, c2, c4] = 1e6
                    count_3D += 1
                    c4 += 1

                elif (fixOMDOT and expOMDOT and np.fabs(omdot_full - omdot) / omdoterr > tolerance):

                    chisq_3D[c1, c2, c4] = 1e6
                    count_3D += 1
                    c4 += 1
                else:

                    if useTEMPO:
    
                        p3 = Popen(com3, stdout=PIPE)
                        out3, err3 = p3.communicate()
                        try:
                            m = re.search(r'Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', out3)
                            chisq_3D[c1,c2,c4] = np.float(m.group(3)) * np.float(m.group(2))
                            if (chisq_3D[c1, c2, c4] < chisq_min):
                                chisq_min = chisq_3D[c1, c2, c4]
                        except:
                            chisq_3D[c1,c2,c4] = 1e6
    
                        c4 += 1
    
                    else:
                        call(com2,stdout=PIPE)
                        p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
                        out3, err3 = p3.communicate()
                        line = out3.split()
                        if (c1 == 0 and c2 == 0 and c4 == 0):
                            nfree = np.float(line[2])
                        chisq_3D[c1,c2, c4] = np.float(line[1]) * np.float(line[2])
                        c4 += 1

                    count_3D += 1

        else:
            fout.close()

            # now run tempo on new/dummy parfile, extract numbers.
            if (useTEMPO):

                if (gridDDGR and (mtot_elem - m2_elem) < 0.):
                    chisq[c1,c2] = 1e6
                    c1 += 1
                    c3 += 1
                    continue

                else:
                     
                    if (fixOMDOT and expOMDOT and np.fabs(omdot - omdot_new) / omdoterr > tolerance):
                        chisq[c1, c2] = 1e6

                    elif (fixGAMMA and expOMDOT and np.fabs(gamma - gamma_new) / gammaerr > tolerance):
                        chisq[c1, c2] = 1e6

                    else: 

                        p = Popen(com2, stdout=PIPE)
                        out, err = p.communicate()
        
                        try:
                            m = re.search('Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', str(out))
                            chisq[c1,c2] = np.float(m.group(3)) * np.float(m.group(2))
                            if (chisq[c1, c2] < chisq_min):
                                chisq_min = chisq[c1, c2]
                        except:
                            chisq[c1,c2] = 1e6
            else:

                if (fixOMDOT and expOMDOT and np.fabs(omdot - omdot_new) / omdoterr > tolerance):
                    chisq[c1, c2] = 1e6

                elif (fixGAMMA and expOMDOT and np.fabs(gamma - gamma_new) / gammaerr > tolerance):
                    chisq[c1, c2] = 1e6

                else:
                    #call(com1, stdout=PIPE)
                    p2 = Popen(com1, stdout=PIPE)
                    out2, err2 = p2.communicate()
                    p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
                    out3, err3 = p3.communicate()
                    line = out3.split()
                    if (c1 == 0 and c2 == 0):
                        nfree = np.float(line[2])
                    chisq[c1,c2] = np.float(line[1]) * np.float(line[2])
                    if (chisq[c1, c2] < chisq_min):
                        chisq_min = chisq[c1, c2]
        c1 += 1
        c3 += 1
    c2 += 1

print()

# create x/y grids for contour plots.
outf = ''

if fixXDOT:
    outf = 'SDgrid.'+obj+'.rs.fixXDOT.pickle'
    
    delta_chisq3D = chisq_3D - np.min(chisq_3D)
    pdf3D = 0.5 * np.exp(-0.5 * delta_chisq3D)

    pdf2D_m2cosi = np.zeros((Ngrid, Ngrid))
    pdf2D_h3stig = np.zeros((Ngrid, Ngrid))
    pdf2D_thetacosi = np.zeros((Ngrid, Ngrid))
    pdf2D_thetastig = np.zeros((Ngrid, Ngrid))

    # loop over two dimensions and collapse to get 2D PDFs.
    for ii in range(Ngrid):
        for jj in range(Ngrid):
            if grid_m2cosi:
                pdf2D_m2cosi[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_thetacosi[jj, ii] = np.sum(pdf3D[:, ii, jj])
            elif gridH3STIG:
                pdf2D_h3stig[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_thetastig[jj, ii] = np.sum(pdf3D[:, ii, jj])

    if grid_m2cosi:
        pdf2D_m2cosi = pdf2D_m2cosi / np.max(pdf2D_m2cosi) * 0.5
        pdf2D_thetacosi = pdf2D_thetacosi / np.max(pdf2D_thetacosi) * 0.5

        plt.figure(1)
        plt.pcolormesh(cosi, m2, pdf2D_m2cosi, vmin=0, vmax=np.max(pdf2D_m2cosi))
        plt.savefig('rs.'+obj+'.THETA.m2cosi.png', format='png')
        plt.figure(2)
        plt.pcolormesh(cosi, theta, pdf2D_thetacosi, vmin=0, vmax=np.max(pdf2D_thetacosi))
        plt.savefig('rs.'+obj+'.THETA.thetacosi.png', format='png')

    elif gridH3STIG:
        pdf2D_h3stig = pdf2D_h3stig / np.max(pdf2D_h3stig) * 0.5
        pdf2D_thetastig = pdf2D_thetastig / np.max(pdf2D_thetastig) * 0.5

        plt.figure(1)
        plt.pcolormesh(stig, h3, pdf2D_h3stig, vmin=0, vmax=np.max(pdf2D_h3stig))
        plt.xlabel(r'$\varsigma$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
        plt.savefig('ortho.'+obj+'.fixTHETA.h3stig.png', format='png')
        plt.figure(2)
        plt.pcolormesh(stig, theta, pdf2D_thetastig, vmin=0, vmax=np.max(pdf2D_thetastig))
        plt.xlabel(r'$\varsigma$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$\Theta$ (deg)', fontproperties=font, fontsize=15)
        plt.savefig('ortho.'+obj+'.fixTHETA.thetastig.png', format='png')

elif fixPX:

    outf = 'SDgrid.'+obj+'.rs.fixPX.pickle'
    
    delta_chisq3D = chisq_3D - np.min(chisq_3D)
    pdf3D = 0.5 * np.exp(-0.5 * delta_chisq3D)

    pdf2D_m1m2 = np.zeros((Ngrid, Ngrid))
    pdf2D_m2cosi = np.zeros((Ngrid, Ngrid))
    pdf2D_h3stig = np.zeros((Ngrid, Ngrid))
    pdf2D_pxm1 = np.zeros((Ngrid, Ngrid))
    pdf2D_pxcosi = np.zeros((Ngrid, Ngrid))
    pdf2D_pxstig = np.zeros((Ngrid, Ngrid))

    # loop over two dimensions and collapse to get 2D PDFs.
    for ii in range(Ngrid):
        for jj in range(Ngrid):
            if grid_m2cosi:
                pdf2D_m2cosi[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_pxcosi[jj, ii] = np.sum(pdf3D[:, ii, jj])
            elif gridM1M2:
                pdf2D_m1m2[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_pxm1[jj, ii] = np.sum(pdf3D[:, ii, jj])
            elif gridH3STIG:
                pdf2D_h3stig[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_pxstig[jj, ii] = np.sum(pdf3D[:, ii, jj])

    if grid_m2cosi:
        pdf2D_m2cosi = pdf2D_m2cosi / np.max(pdf2D_m2cosi) * 0.5
        pdf2D_pxcosi = pdf2D_pxcosi / np.max(pdf2D_pxcosi) * 0.5

        plt.figure(1)
        plt.pcolormesh(cosi, m2, pdf2D_m2cosi, vmin=0, vmax=np.max(pdf2D_m2cosi))
        plt.savefig('rs.'+obj+'.fixPX.m2cosi.png', format='png')
        plt.xlabel(r'$\cos i$', fontproperties=font, fontsize=15)
        plt.ylabel(r'Companion Mass (M$_\odot$)', fontproperties=font, fontsize=15)
        plt.figure(2)
        plt.pcolormesh(cosi, px, pdf2D_pxcosi, vmin=0, vmax=np.max(pdf2D_pxcosi))
        plt.xlabel(r'$\cos i$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$\varpi$ (mas)', fontproperties=font, fontsize=15)
        plt.savefig('rs.'+obj+'.fixPX.pxcosi.png', format='png')

    elif gridM1M2:
        pdf2D_m1m2 = pdf2D_m1m2 / np.max(pdf2D_m1m2) * 0.5
        pdf2D_pxm1 = pdf2D_pxm1 / np.max(pdf2D_pxm1) * 0.5

        plt.figure(1)
        plt.pcolormesh(m1, m2, pdf2D_m1m2, vmin=0, vmax=np.max(pdf2D_m1m2))
        plt.savefig('rs.'+obj+'.fixPX.m1m2.png', format='png')
        plt.xlabel(r'Pulsar Mass (M$_\odot$)', fontproperties=font, fontsize=15)
        plt.ylabel(r'Companion Mass (M$_\odot$)', fontproperties=font, fontsize=15)
        plt.figure(2)
        plt.pcolormesh(m1, px, pdf2D_pxm1, vmin=0, vmax=np.max(pdf2D_pxm1))
        plt.xlabel(r'Pulsar Mass (M$_\odot$)', fontproperties=font, fontsize=15)
        plt.ylabel(r'$\varpi$ (mas)', fontproperties=font, fontsize=15)
        plt.savefig('rs.'+obj+'.fixPX.pxm1.png', format='png')


    elif gridH3STIG:
        pdf2D_h3stig = pdf2D_h3stig / np.max(pdf2D_h3stig) * 0.5
        pdf2D_pxstig = pdf2D_pxstig / np.max(pdf2D_pxstig) * 0.5

        plt.figure(1)
        plt.pcolormesh(stig, h3, pdf2D_h3stig, vmin=0, vmax=np.max(pdf2D_h3stig))
        plt.xlabel(r'$\varsigma$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
        plt.savefig('ortho.'+obj+'.fixPX.h3stig.png', format='png')
        plt.figure(2)
        plt.pcolormesh(stig, px, pdf2D_pxstig, vmin=0, vmax=np.max(pdf2D_pxstig))
        plt.xlabel(r'$\varsigma$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$\varpi$ (mas)', fontproperties=font, fontsize=15)
        plt.savefig('ortho.'+obj+'.fixPX.pxstig.png', format='png')


else:
    deltachi2 = chisq - np.min(chisq)
    pdf2D = 0.5 * np.exp(-0.5 * deltachi2)

    outf_contour = ''

    if (fixOMDOT):
        outf_contour = 'rs.'+obj+'.fixOMDOT.png'
    else:
        outf_contour = 'rs.'+obj+'.png'

    plt.pcolormesh(x, y, pdf2D, vmin=0, vmax=np.max(pdf2D), cmap="Blues", shading="auto")
    plt.colorbar()

    if gridH3STIG:
        outf_contour = 'h3stig.'+obj+'.orig.png'
        plt.xlabel(r'$\varsigma$', fontproperties=font, fontsize=15)
        plt.ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
    elif gridH3H4:
        outf_contour = 'h3h4.'+obj+'.orig.png'
        plt.xlabel(r'$h_4$ ($\mu$s)', fontproperties=font, fontsize=15)
        plt.ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
    elif gridM1M2:
        plt.xlabel(r'Pulsar Mass (${\rm M}_\odot$)', fontproperties=font, fontsize=15)
        plt.ylabel(r'Companion Mass (M$_\odot$)', fontproperties=font, fontsize=15)
    elif gridM2MTOT:
        plt.xlabel(r'Total Mass (${\rm M}_\odot$)', fontproperties=font, fontsize=15)
        plt.ylabel(r'Companion Mass (M$_\odot$)', fontproperties=font, fontsize=15)
    else:
        plt.xlabel(r'$\cos i$', fontproperties=font, fontsize=15)
        plt.ylabel(r'Companion Mass (M$_\odot$)', fontproperties=font, fontsize=15)
    plt.savefig(outf_contour, format='png')

    outf = ''

    if (fixOMDOT):
        outf = 'SDgrid.'+obj+'.rs.fixOMDOT.pickle'
    else:
        outf = 'SDgrid.'+obj+'.rs.pickle'

GridDict = {}
GridDict['PSR'] = obj

if grid_m2cosi:
    GridDict['COSI'] = cosi
    GridDict['M2']   = m2

elif gridH3H4:
    GridDict['Tsun'] = T_sun
    GridDict['H4'] = h4
    GridDict['H3'] = h3

elif gridH3STIG:
    GridDict['Tsun'] = T_sun
    GridDict['STIG'] = stig
    GridDict['H3']   = h3

elif (gridDDGR or gridM2MTOT):
    GridDict['MTOT'] = mtot
    GridDict['M2'] = m2

else:
    GridDict['M1'] = m1
    GridDict['M2'] = m2
   

GridDict['massfunc'] = massfunc

if fixXDOT:

    GridDict['THETA'] = theta
    GridDict['chisq'] = chisq_3D

elif fixPX:

    GridDict['PX'] = px
    GridDict['PX_bestfit'] = px_bestfit
    GridDict['chisq'] = chisq_3D

else:

    GridDict['chisq'] = chisq

print(np.min(chisq), np.max(chisq))

GridDict['chisq_bestfit'] = chisq_bestfit
GridDict['M2_bestfit'] = m2_bestfit
GridDict['SINI_bestfit'] = sini_bestfit
GridDict['DOF'] = nfree

# save variables for future use.
pout = open(outf, 'wb')
pickle.dump(GridDict, pout)
pout.close()

# clean up shop.
q = Popen(['rm','new.par','dummygrid.par'],stdout=PIPE)
