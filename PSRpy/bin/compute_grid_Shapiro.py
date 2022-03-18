#! /usr/bin/python

from astropy.coordinates import SkyCoord
from subprocess import call, Popen, PIPE
from PSRpy.parfile import Parfile
from astropy import units as u
from PSRpy.const import T_sun
from os.path import isfile
from copy import deepcopy
import PSRpy.orbit.variations as orbvar
import PSRpy.orbit.masses as masses
import PSRpy.orbit.pkcorr as pkcorr
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
    "-G", action="store_true", dest="useGLS", 
    help="Use GLS fitting in TEMPO."
)

parser.add_argument(
    "-n", action="store", dest="Ngrid", default=50, type=int, 
    help="Number of grid points to generate. (Default: 50)"
)

parser.add_argument(
    "-p", action="store", dest="pnum", default=None, type=str, 
    help="TEMPO pulse-number file."
)

parser.add_argument(
    "--M2", nargs=2, action="store", metavar=("min", "max"), 
    dest="M2limits", default=[0.01, 1.], 
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
    "--M1M2", action="store_true", dest="grid_m1m2", 
    help='Uniformly grid over M1/M2.'
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

# decide which type of grid to generate.
grid_m2cosi = True
x = None
y = None
z = None
x_label = ""
y_label = ""
z_label = ""

if any([gridDDGR, gridDDK, gridM2MTOT, gridM1M2, gridH3STIG, gridH3H4]):
    grid_m2cosi = False

    if any([gridDDGR, gridM2MTOT]):
        x_label = "MTOT"
        y_label = "M2"
        x = np.linspace(mtot_lo, mtot_hi, num=Ngrid)
        y = np.linspace(m2_lo, m2_hi, num=Ngrid)
    
    elif gridDDK:
        x_label = r"\cos i"
        y_label = "KOM"
        x = np.linspace(cosi_lo, cosi_hi, num=Ngrid)
        y = np.linspace(theta_lo, theta_hi, num=Ngrid)

    elif gridM1M2:
        x_label = "M1"
        y_label = "M2"
        x = np.linspace(m1_lo, m1_hi, num=Ngrid)
        y = np.linspace(m2_lo, m2_hi, num=Ngrid)

    elif gridH3STIG:
        x_label = "STIG"
        y_label = "H3"
        x = np.linspace(stig_lo, stig_hi, num=Ngrid)
        y = np.linspace(h3_lo, h3_hi, num=Ngrid)

    elif gridH3H4:
        x_label = "H4"
        y_label = "H3"
        x = np.linspace(h4_lo, h4_hi, num=Ngrid)
        y = np.linspace(h3_lo, h3_hi, num=Ngrid)
        
    else:
        sys.exit("ERROR: can't recognize grid and determine axes!")

else:
    x_label = "COSI"
    y_label = "M2"
    x = np.linspace(cosi_lo, cosi_hi, num=Ngrid)
    y = np.linspace(m2_lo, m2_hi, num=Ngrid)


# determine if desired grid is over three coordinates.
grid_3D = False

if any([px_lo, px_hi]) != 0.:
    fixPX = True
    grid_3D = True
    z_label = "PX"
    z = np.linspace(px_lo, px_hi, num=Ngrid)

elif grid_m2cosi and any([theta_lo, theta_hi]) != 0.:
    fixXDOT = True
    grid_3D = True
    z_label = "THETA"
    z = np.linspace(theta_lo, theta_hi, num=Ngrid)

elif any([xomdot_lo, xomdot_hi]) != 0.:
    fixXOMDOT = True
    grid_3D = True
    z_label = "XOMDOT"
    z = np.linspace(xomdot_lo, xomdot_hi, num=Ngrid)

### set grid  parameters.
chisq = None 

if grid_3D:
    chisq = np.zeros((Ngrid, Ngrid, Ngrid))

else:
    chisq = np.zeros((Ngrid, Ngrid))

### read in supplied parfile into PSRpy parfile object and extract info.
input_par = Parfile()
input_par.read(inpar)
obj = ""
num_degrees_of_freedom = 0

if input_par.PSR["value"] is not None:
    obj = input_par.PSR["value"]

elif input_par.PSRJ["value"] is not None:
    obj = input_par.PSRJ["value"]

else:
    sys.exit("ERROR: parfile must have PSR or PSRJ set!")

# extract Keplerian elements, compute mass function.
eccentricity = input_par.E["value"]
orbital_period = input_par.PB["value"]
semimajor_axis = input_par.A1["value"]
mass_func = masses.mass_function(orbital_period, semimajor_axis)

if eccentricity is None:
    if input_par.EPS1["value"] is not None and input_par.EPS2["value"] is not None:
        eccentricity = np.sqrt(input_par.EPS1["value"]**2 + input_par.EPS2["value"]**2)
 
    else:
        sys.exit("ERROR: input parfile seems to have incomplete binary model?")

# extract position info and compute converted Galactic coordinates.
longitude = None
latitude = None
units = None
frame = ""

if input_par.RAJ["value"] is not None and input_par.DECJ["value"] is not None:
    latitude = input_par.DECJ["value"]
    longitude = input_par.RAJ["value"]
    units = (u.hourangle, u.deg)
    frame = "icrs"

elif input_par.LAMBDA["value"] is not None and input_par.BETA["value"] is not None:
    latitude = input_par.BETA["value"]
    longitude = input_par.LAMBDA["value"]
    units = (u.deg, u.deg)
    frame = "barycentrictrueecliptic"

else:
    sys.exit("ERROR: astrometric coordinates are missing or inconsistent!")

coords = SkyCoord(longitude, latitude, frame=frame, unit=units)

# if available, compute total magnitude of proper motion.
pm = 0.

if input_par.PMRA["value"] is not None and input_par.PMDEC["value"] is not None:
    pm = np.sqrt(input_par.PMRA["value"]**2 + input_par.PMDEC["value"]**2)

elif input_par.PMLAMBDA["value"] is not None and input_par.PMBETA["value"] is not None:
    pm = np.sqrt(input_par.PMLAMBDA["value"]**2 + input_par.PMBETA["value"]**2)

else:
    print("WARNING: both proper motion terms are not in parfile!")
    print("... neglecting all calculations that involve these terms ...")

# before proceeding, convert proper motion to rad/s.
pm *= (np.pi / 180 / 1000 / 3600 / 365.25 / 86400)

### determine structure of form for TEMPO or TEMPO2 call, depending
### on inputs supplied at the command line.
com1, com2 = [], []
chisq_bestfit = 0.
chisq_min = 1e6

# if TEMPO is desired, determine the call.
if useTEMPO:

    # if GLS is desired, add the "-G" flag.
    if useGLS:

        # if pulse numbers are to be used, add the flag and pulse-number file.
        if inpnum is not None:
            com1 = ['tempo', '-ni', inpnum, '-G', '-f', inpar, intim]
            com2 = ['tempo', '-ni', inpnum, '-G', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-ni', inpnum, '-G', '-f', './dummygrid2.par', intim]
        
        # otherwise, ignore pulse numbers.
        else:
            com1 = ['tempo', '-G', '-f', inpar, intim]
            com2 = ['tempo', '-G', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-G', '-f', './dummygrid2.par', intim]

    # same as above if-statement, but neglecting GLS fit option.
    else:

        if inpnum is not None:
            com1 = ['tempo', '-ni', inpnum, '-f', inpar, intim]
            com2 = ['tempo', '-ni', inpnum, '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-ni', inpnum, '-f', './dummygrid2.par', intim]

        else:
            com1 = ['tempo', '-f', inpar, intim]
            com2 = ['tempo', '-f', './dummygrid.par', intim]
            com3 = ['tempo', '-f', './dummygrid2.par', intim]

# if TEMPO2 is instead desired, suss out the call.
else:

    # if QR composition is desired, add the right flag.
    if useQRFIT:
        com1 = ['tempo2', '-qrfit', '-newpar', '-f', './dummygrid.par', intim]
        com2 = ['tempo2', '-qrfit', '-newpar', '-f', './dummygrid2.par', intim]

    else:
        com1 = ['tempo2', '-newpar', '-f', './dummygrid.par', intim]
        com2 = ['tempo2', '-newpar', '-f', './dummygrid2.par', intim]

# before performing computations, print relevant info to terminal.
print("Getting Shapiro-delay grid...")

if useGLS:
    print("    * using GLS...")

if fixOMDOT:
    print("    * fixing OMDOT to GR value", end="")

    # check if input parfile has an OMDOT.
    if input_par.OMDOT["value"] is not None:
        print("(best-fit OMDOT = {0:.5f} deg/yr)".format(input_par.OMDOT["value"]))

    else:
        print()

if fixPBDOT:

    # check if input parfile has an PBDOT.
    if input_par.PBDOT["value"] is not None and fixPX:
        print("    * fixing PBDOT to GR+kinematic value (best-fit PBDOT = {0:.5f} * 1e-12)".format(
                input_par.PBDOT["value"]
            )
        )

    elif hasattr(input_par, "PBDOT"):
        print("    * fixing PBDOT to GR value (best-fit PBDOT = {0:.5f} * 1e-12)".format(
                input_par.PBDOT
            )
        )

    else:
        sys.exit("ERROR: no PBDOT in original parfile; must add an PBDOT line in parfile!")

if fixGAMMA:

    # check if input parfile has a GAMMA.
    if hasattr(input_par, "GAMMA"):
        print("    * fixing GAMMA to GR value (best-fit GAMMA = {0:.5f} ms)".format(
                input_par.GAMMA["value"] * 1000.
            )
        )

    else:
        sys.exit("ERROR: no GAMMA in original parfile; must add an GAMMA line in parfile!")

if fixXDOT:
    print("Fixing XDOT (best-fit XDOT = {0:.5f}) * 1e-12".format(input_par.XDOT["value"]))
    cosi_bestfit = np.sqrt(1 - input_par.SINI["value"]**2)
    coti_bestfit = cosi_bestfit / input_par.SINI["value"]
    xdot_max = pm * semimajor_axis * coti_bestfit * 1e12
    sinTHETA = input_par.XDOT["value"] / xdot_max
    print("    * predicted max. XDOT = {0:.5f}...".format(xdot_max))
    print("    * predicted THETA     = {0:.5f}...".format(np.arcsin(sinTHETA) * 180 / np.pi % 360))

# if desired, print info on orthometric parameters.
if gridH3H4:
    print("    * orthometric (approximate) timing solution...")
    print("        - best-fit H3: {0:.3f} microsec".format(input_par.H3 * 1e6))
    print("        - best-fit H4: {0:.3f} microsec".format(input_par.H4 * 1e6))

if gridH3STIG:
    print("    * orthometric (full) timing solution...")
    print("        - best-fit H3:   {0:.3f} microsec".format(input_par.H3 * 1e6))

    # the inclination parameter is called different names in TEMPO and TEMPO2, 
    # so treat this calculation carefully.
    stig_bestfit = 0.

    if hasattr(input_par, "STIG"):
        stig_bestfit = input_par.STIG

    elif hasattr(input_par, "VARSIGMA"):
        stig_bestfit = input_par.VARSIGM

    print("        - best-fit STIG: {0:.3f}".format(stig_bestfit))
    sini_bestfit = 2 * stig_bestfit / (1 + stig_bestfit**2)
    m2_bestfit = input_par.H3 / stig_bestfit**3 / T_sun.value
    print("    * derived M2:    {0:.3f} solar masses".format(m2_bestfit))
    print("    * derived SINI:  {0:.3f}".format(sini_bestfit))

if gridDDGR:
    print("Fixing MTOT (best-fit = {0:.5f})...".format(input_par.MTOT["value"]))

# define variables to be re-defined at each grid point if applicable
# to the desired grid type.
c2, c3 = 0, 0
new_par = deepcopy(input_par)
gamma_new = 0.
pbdot_elem = 0.
pbdot_full = 0.
xdot_new = 0.
omdot_elem = 0.
omdot_full = 0.
sini_elem = 0.
cosi_elem = 0.
shapmax_elem = 0.
m1_elem = 0.
m2_elem = 0.
mtot_elem = 0.
count_3D = 0

### now perform loops over the axes for grid calculation.
# first, loop over x coordinate.
for x_elem in x:

    # define appropriate variables, depending on desired grid.
    if grid_m2cosi:
        cosi_elem = x_elem
        sini_elem = np.sqrt(1 - cosi_elem**2)
        shapmax_elem = -np.log(1 - sini_elem)

    elif gridDDGR or gridM2MTOT:
        mtot_elem = x_elem

    elif gridM1M2:
        m1_elem = x_elem

    c1 = 0

    # now loop over y coordinate.
    for y_elem in y:

        # if grid uses "traditional" parameters, compute M1/M2/MTOT/SINI values.
        if any([grid_m2cosi, gridDDGR, gridM1M2, gridM2MTOT]):
            m2_elem = y_elem
            new_par.set("M2", {"value": m2_elem, "flag": 0})

            # if other parameter is not MTOT, then set the correct parameter.
            if grid_m2cosi:
                m1_elem = np.sqrt((m2_elem * sini_elem)**3 / mass_func) - m2_elem
                mtot_elem = m1_elem + m2_elem

                if gridDDK:
                    new_par.set("KIN", {"value": np.arccos(cosi_elem) * 180 / np.pi, "flag": 0})

                else:
                    new_par.set("SINI", {"value": sini_elem, "flag": 0})

            # otherwise, set MTOT and define other parameters related to it.
            elif any([gridDDGR, gridM1M2, gridM2MTOT]):
                new_par.set("MTOT", {"value": mtot_elem, "flag": 0})

                if (gridDDGR or gridM2MTOT):
                    m1_elem = mtot_elem - m2_elem
                
                else: 
                    mtot_elem = m1_elem + m2_elem

                sini_elem = (mass_func * (m1_elem + m2_elem)**2)**(1./3.) / m2_elem

        # if grid is instead over KOM/KIN terms, treat those separately.
        elif gridDDK:
            new_par.set("KIN", {"value": np.arccos(x_elem) * 180 / np.pi, "flag": 0})
            new_par.set("KOM", {"value": y_elem, "flag": 0})

        # if grid is instead "orthometric", compute H3/H4/STIG values.
        elif any([gridH3H4, gridH3STIG]):
            pass

        # if desired, compute corresponding post-Keplerian parameters here.
        # the following is for OMDOT.
        if any([fixOMDOT, fixXOMDOT, fixXDOT]):
            omdot_elem = orbvar.omdot_GR(
                m1_elem, m2_elem, orbital_period, eccentricity, use_PK2=use_PK2
            )

            if fixOMDOT:
                new_par.set("OMDOT", {"value": omdot_elem, "flag": 0})

        # the following is for GAMMA.
        if fixGAMMA:
            gamma_elem = orbvar.gamma_GR(
                m1_elem, m2_elem, orbital_period, eccentricity
	    )

            new_par.set("GAMMA", {"value": gamma_elem, "flag": 0})

        # the following is for PBDOT.
        if fixPBDOT:
            pbdot_elem = orbvar.pbdot_GR(
                m1_elem, m2_elem, orbital_period, eccentricity
            )

            # if not gridding over PX/DIST, then just compute GR term and set.
            if not fixPX:
                new_par.set("PBDOT", {"value": pbdot_elem, "flag": 0})

        # if grid is 3D, do not yet update terminal info on percent completed.
        if fixXDOT or fixPX:
            pass

        # otherwise, print percent completed to terminal screen for a 2D grid.
        else:
            perc = 100. * (c3 + 1.) / Ngrid**2
            sys.stdout.write('\rPercent completed: {0:.3f}% (min. chisq = {1:.3f})'.format(
                    perc, chisq_min
                )
            )
            sys.stdout.flush()

        # if grid is 3D, loop over third dimension.
        # first consider if grid over Kopeikin terms is desired.
        if z is not None:
            c4 = 0

            for z_elem in z:
                perc = 100. * (count_3D + 1) / Ngrid**3
                sys.stdout.write(
                        "\rPercent completed: {0:.3f}% (min. chisq = {1:.3f})".format(
                        perc, chisq_min
                    )
                )
                sys.stdout.flush()
                c3 += 1

                if fixXDOT and not gridDDK:
                    coti_elem = np.sqrt(1 - sini_elem**2) / sini_elem
                    xdot_full = semimajor_axis * pm * coti_elem * np.sin(z_elem * np.pi / 180)  
                    omdot_sec = pm / sini_elem * np.cos(z_elem * np.pi / 180)
                    omdot_sec *= (180 / np.pi * 86400 * 365.25)
                    omdot_full = omdot_elem + omdot_sec

                    # write out XDOT in correct units if using TEMPO or TEMPO2.
                    if useTEMPO:
                        new_par.set("XDOT", {"value": xdot_full * 1e12, "flag": 0})
                    else:
                        new_par.set("XDOT", {"value": xdot_full, "flag": 0})

                    # if model is not ELL1, write out OMDOT.
                    if (new_par.BINARY["value"] != "ELL1"):
                        new_par.set("OMDOT", {"value": omdot_full, "flag": 0})

                # else, assume DDK model and set KOM.
                elif gridDDK:
                    pass
                    #fout2.write("{0}                 {1:.8f}   0\n".format('KOM', theta_elem))

                # if grid is over PX/DIST, compute relevant terms and set.
                elif fixPX:
                    new_par.set("PX", {"value": z_elem, "flag": 0})
                    
                    if fixPBDOT:
                        dop_comps, err = pkcorr.doppler(
                            1 / z_elem, 0.1, coords.galactic.b.deg, coords.galactic.l.deg, 
                            pm, 0.1
                        )
                        pbdot_full = pbdot_elem + np.sum(dop_comps) * orbital_period * 86400
                        new_par.set("PBDOT", {"value": pbdot_full * 1e12, "flag": 0})

                # if grid is over XOMDOT, compute relevant terms and set.
                elif fixXOMDOT:
                    omdot_full = omdot_elem + z_elem * 1e-6
                    new_par.set("OMDOT", {"value": omdot_full, "flag": 0})

                # write temporary parfile to disk.
                new_par.write(outfile="dummygrid2.par")

                # finally, determine best-fit chisq at current grid point.
                if useTEMPO:
                    p3 = Popen(com3, stdout=PIPE)
                    out3, err3 = p3.communicate()

                    try:
                        m = re.search(
                            r'Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', str(out3)
                        )
                        chisq[c1,c2,c4] = float(m.group(3)) * int(m.group(2))

                        if (chisq[c1, c2, c4] < chisq_min):
                            chisq_min = chisq[c1, c2, c4]
                            num_degrees_of_freedom = int(m.group(2))

                    except:
                        chisq[c1,c2,c4] = 1e6

                    c4 += 1

                else:
                    call(com2,stdout=PIPE)
                    p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
                    out3, err3 = p3.communicate()
                    line = out3.split()
                    if (c1 == 0 and c2 == 0 and c4 == 0):
                        num_degrees_of_freedom = float(line[2])
                    chisq[c1,c2, c4] = np.float(line[1]) * np.float(line[2])
                    c4 += 1

                count_3D += 1

        # if instead a 2D grid is desired, just proceed to calling tempo/tempo2.
        else:
            new_par.write(outfile="dummygrid.par")

            # if desired, use tempo for chisq estimation.
            if useTEMPO:
                p = Popen(com2, stdout=PIPE)
                out, err = p.communicate()

                try:
                    m = re.search('Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', str(out))
                    chisq[c1, c2] = float(m.group(3)) * float(m.group(2))
                    nfree = int(m.group(2))

                    # if current chisq is lowest yet, record.
                    if (chisq[c1, c2] < chisq_min):
                        chisq_min = chisq[c1, c2]

                except:
                    chisq[c1,c2] = 1e6

            else:
                #call(com1, stdout=PIPE)
                p2 = Popen(com1, stdout=PIPE)
                out2, err2 = p2.communicate()
                p3 = Popen(['grep','CHI2R','./new.par'],stdout=PIPE)
                out3, err3 = p3.communicate()
                line = out3.split()
                if (c1 == 0 and c2 == 0):
                    nfree = np.float(line[2])
                chisq[c1,c2] = float(line[1]) * float(line[2])
                if (chisq[c1, c2] < chisq_min):
                    chisq_min = chisq[c1, c2]

        # update indices for y coordinate, bookkeeping.
        c1 += 1
        c3 += 1

    # update index for x coordinate.
    c2 += 1

print()
print(np.min(chisq), np.max(chisq))

# before writing data to disk, create a plot or two.
outf = ''

if z is not None:
    delta_chisq3D = chisq - np.min(chisq)
    pdf3D = 0.5 * np.exp(-0.5 * delta_chisq3D)
    pdf2D_xy = np.zeros((Ngrid, Ngrid))
    pdf2D_xz = np.zeros((Ngrid, Ngrid))
    pdf2D_yz = np.zeros((Ngrid, Ngrid))
    outf_contour_xy = ""
    outf_contour_xz = ""
    outf_contour_yz = ""

    if fixOMDOT:
        outf_contour_xy = "rs.{0}.{1}{2}.fixOMDOT.png".format(obj, x_label, y_label)
        outf_contour_xz = "rs.{0}.{1}{2}.fixOMDOT.png".format(obj, x_label, z_label)
        outf_contour_yz = "rs.{0}.{1}{2}.fixOMDOT.png".format(obj, y_label, z_label)
    else:
        outf_contour_xy = "rs.{0}.{1}{2}.png".format(obj, x_label, y_label)
        outf_contour_xz = "rs.{0}.{1}{2}.png".format(obj, x_label, z_label)
        outf_contour_yz = "rs.{0}.{1}{2}.png".format(obj, y_label, z_label)

    # loop over two dimensions and collapse to get 2D PDFs.
    for ii in range(Ngrid):
        for jj in range(Ngrid):
                pdf2D_xy[jj, ii] = np.sum(pdf3D[jj, ii, :])
                pdf2D_xz[jj, ii] = np.sum(pdf3D[jj, :, ii])
                pdf2D_yz[jj, ii] = np.sum(pdf3D[:, jj, ii])

    plt.pcolormesh(x, y, pdf2D_xy, vmin=0, vmax=np.max(pdf2D_xy), cmap="Blues", shading="auto")
    plt.colorbar()
    plt.xlabel(x_label, fontproperties=font, fontsize=15)
    plt.ylabel(y_label, fontproperties=font, fontsize=15)
    plt.savefig(outf_contour_xy, fmt="png")
    plt.clf()

    plt.pcolormesh(x, z, pdf2D_xz, vmin=0, vmax=np.max(pdf2D_xz), cmap="Blues", shading="auto")
    plt.colorbar()
    plt.xlabel(x_label, fontproperties=font, fontsize=15)
    plt.ylabel(z_label, fontproperties=font, fontsize=15)
    plt.savefig(outf_contour_xz, fmt="png")
    plt.clf()

    plt.pcolormesh(y, z, pdf2D_yz, vmin=0, vmax=np.max(pdf2D_yz), cmap="Blues", shading="auto")
    plt.colorbar()
    plt.xlabel(y_label, fontproperties=font, fontsize=15)
    plt.ylabel(z_label, fontproperties=font, fontsize=15)
    plt.savefig(outf_contour_yz, fmt="png")

else:
    deltachi2 = chisq - np.min(chisq)
    pdf2D = 0.5 * np.exp(-0.5 * deltachi2)
    outf_contour = ""

    if fixOMDOT:
        outf_contour = "rs.{0}.fixOMDOT.png".format(obj)
    else:
        outf_contour = "rs.{0}.png".format(obj)

    plt.pcolormesh(x, y, pdf2D, vmin=0, vmax=np.max(pdf2D), cmap="Blues", shading="auto")
    plt.colorbar()
    plt.xlabel(x_label, fontproperties=font, fontsize=15)
    plt.ylabel(y_label, fontproperties=font, fontsize=15)
    plt.savefig(outf_contour, fmt="png")

### finally, load data dictionary and write to file.
GridDict = {}
GridDict['PSR'] = obj
GridDict['Tsun'] = T_sun
GridDict['massfunc'] = mass_func
GridDict['chisq'] = chisq
GridDict['chisq_bestfit'] = np.min(chisq)
GridDict[x_label] = x
GridDict[y_label] = y

if z is not None:
    GridDict[z_label] = z

GridDict['M2_bestfit'] = input_par.M2["value"]
GridDict['PX_bestfit'] = input_par.PX["value"]
GridDict['DOF'] = num_degrees_of_freedom

# be careful about how to stash the best-fit value of SINI
if input_par.KIN["value"] is not None:
    GridDict['SINI_bestfit'] = np.sin(input_par.KIN["value"] * 180 / np.pi)

else:
    GridDict['SINI_bestfit'] = input_par.SINI["value"]

# save variables for future use.
pout = open("SDgrid.{0}.pkl".format(obj), 'wb')
pickle.dump(GridDict, pout)
pout.close()

# clean up shop.
q = Popen(['rm','new.par','dummygrid.par'],stdout=PIPE)
