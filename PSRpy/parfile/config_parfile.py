"""
Configuration file containing parameters and lists that control the Parfile() object.
"""

# variables that control lengths of parameter lists.
n_derivatives_orbit = 10
n_derivatives_spin = 20
n_derivatives_DM = 10
n_derivatives_FD = 6
n_orbits_BTX = 3
n_bins_DMX = 500

# list of astrometric parameters.
parameter_list_astrometry = ["RAJ", "DECJ", "PMRA", "PMDEC"] +\
    ["LAMBDA", "BETA", "PMLAMBDA", "PMBETA"] +\
    ["ELON", "ELAT", "PMELON", "PMELAT"] +\
    ["PX"]

# list of spin-frequency derivatives.
parameter_list_spin = ["F{0}".format(str(num)) for num in range(n_derivatives_spin)]

# lst of DM parameters.
parameter_list_DM = ["DM"] +\
    ["DM{0}".format(num) for num in range(1, n_derivatives_DM)]

parameter_list_DMX = ["DMX"]
for current_bin in range(n_bins_DMX):
    for current_parameter in ["DMX", "DMXEP", "DMXR1", "DMXR2", "DMXF1", "DMXF2"]:
        parameter_list_DMX += ["{0}_{1}".format(current_parameter, str(current_bin).zfill(4))]

# list of FD parameters.
parameter_list_FD = ["FD{0}".format(num) for num in range(1, n_derivatives_FD)]

# lists for orbital parameters.
parameter_list_orbit_Kepler = [
    "A1", "PB", "E", "ECC", "OM", "TO", "EPS1", "EPS2", "TASC"
]

parameter_list_orbit_PK = ["M2", "SINI", "GAMMA", "DTHETA", "DR"]

parameter_list_orbit_derivatives = {
    "A1" : ["XDOT"] + ["X{0}".format(str(num) for num in range(2, n_derivatives_orbit))],
    "E"  : ["EDOT"] + ["E{0}".format(str(num) for num in range(2, n_derivatives_orbit))],
    "OM" : ["OMDOT"] + ["OM{0}".format(str(num) for num in range(2, n_derivatives_orbit))],
    "FB" : ["FB{0}".format(str(num) for num in range(1, n_derivatives_orbit))],
    "PB" : ["PBDOT"] # to my knowledge, only PBDOT is defined in TEMPO/TEMPO2/PINT.
}

parameter_list_orbit_BTX = []
for current_orbit in range(2, n_orbits_BTX):
    for current_parameter in ["A1", "PB", "E", "OM", "T0"]:
        parameter_list_orbit_BTX += ["{0}_{1}".format(current_parameter, current_orbit)]

# list of parameters with string valuesl.
parameter_list_string = ["PSR", "PSRJ", "RAJ", "DECJ", "SOLARN0", "EPHEM", "ECL", 
    "CLK", "UNITS", "TIMEEPH", "T2CMETHOD", "CORRECT_TROPOSPHERE", "PLANET_SHAPIRO",
    "DILATEFREQ", "NTOA", "TRES", "TZRMJD", "TZRFRQ", "TZRSITE", "NITS", "INFO",
    "BINARY", "DCOVFILE"
]

# list of non-physical parameters with integer values.
parameter_list_int = ["NTOA", "NITS", "NDDM", "EPHVER", "MODE"]

# list of red/white-noise parameters for TOA and DM uncertainties.
parameter_list_error = ["T2EFAC", "T2EQUAD", "ECORR", "JUMP", "DMEFAC", "DMJUMP"]

# list of all parameters defined above.
parameter_list_full = parameter_list_string + parameter_list_int + parameter_list_error +\
    parameter_list_astrometry +\
    parameter_list_spin +\
    parameter_list_DM + parameter_list_DMX +\
    parameter_list_FD +\
    parameter_list_orbit_Kepler +\
    parameter_list_orbit_BTX

# list of binary models.
model_list_binary_eccentric = ["DD", "DDGR", "DDFWHE", "BT", "BTX"]
model_list_binary_circular = ["ELL1", "ELL1H"]
