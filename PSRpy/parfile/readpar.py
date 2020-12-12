#! /usr/bin/python

from re import match
from PSRpy.const import c, G, M_sun, T_sun
from astropy.coordinates import Angle
from .config_parfile import *
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import sys

def represents_an_int(input_string):
    """
    Checks if string can be represented as a Python integer and returns a boolean.
    This function is used to evaluate if 

    taken from: https://stackoverflow.com/questions/1265665/how-can-i-check-if-a-string-represents-an-int-without-using-try-except
    """

    is_an_integer = False

    try:
        value = int(input_string)
        is_an_interger = True

    except ValueError:
        pass

    return is_an_integer

class Parfile(object):
    """
    Reads in parfile parameter values, uncertainties and flags as object attributes.

    Inputs
    ------

    input_parfile : str
        pulsar-timing parameter file

    efac : float
        constant factor to multiply across all uncertainties.
    
    Notes
    -----
    Parameters stored as class attibutes using same parameter-file name (e.g. par.RAJ 
    stores the 'RAJ' value in 'par' object). Parameter uncertainties are stored as attributes 
    with parameter name plus an 'err' extension (e.g. par.RAJerr stores the 'RAJ' uncertainty 
    in the 'par' object).
    """


    def __init__(self, input_parfile, efac=1):

        # first, read input data.
        self._read_parfile(input_parfile, efac=efac)
        
        # next, derive parameters based on parfile data.
        # TODO: define internal-use function for this purpose.


    def _read_parfile(self, input_parfile, efac=1):
        # preserve order of parameters in parfile.
        self.fit_parameters = []
        parorder = []

        for line in open(input_parfile, "r").readlines():
            lsplit = line.split()

            if (len(lsplit) != 0):
                # if this line is commented out, then skip.
                if (lsplit[0] == 'C' or lsplit[0] == '#'):
                    continue

                # else, loop will continue.
                parname, parvalue = lsplit[0], lsplit[1]

                # set the following attributes as strings.
                if (parname in parameter_list_string):
                    setattr(self, parname, parvalue)
                    # the following is for 'RAJ', 'DECJ' that have flags/errors.
                    if (len(lsplit) > 2):
                        # check if third element is an integer;
                        # if not, it is most likely the parameter uncertainty.
                        # (some ATNF parfiles have this circumstance.)
                        is_flag = represents_an_int(lsplit[2])

                        if (is_flag):
                            setattr(self, parname + 'flag', np.int(lsplit[2]))
                            if (np.int(lsplit[2]) == 1):
                                self.fit_parameters.append(parname)

                        else:
                            setattr(self, parname + 'err', np.float(lsplit[2]))

                    if (len(lsplit) > 3):
                        setattr(self, parname + 'err', efac * np.float(lsplit[3]))

                # set these as integers.
                elif (parname in parameter_list_int):
                    setattr(self, parname, np.int(parvalue))

                # otherwise, if not a JUMP, assume it's a fit parameter 
                # and store value/errors as floats.
                elif (parname not in parameter_list_error):
                    # switch 'D' with 'e' for exponents.
                    if (parvalue.find('D') != -1):
                        parvalue = parvalue.replace('D','e')
                    # DDK model in TEMPO can have SINI = KIN.
                    if (parvalue == 'KIN'):
                        setattr(self, parname, parvalue)
                    else:
                        setattr(self, parname, np.longdouble(parvalue))
                    # set flag and error attributes if present in parfile.
                    if (len(lsplit) > 2):
                        # check if third element is an integer;
                        # if not, it is most likely the parameter uncertainty.
                        # (some ATNF parfiles have this circumstance.)
                        is_flag = represents_an_int(lsplit[2])

                        if (parname != 'START' and parname != 'FINISH' and is_flag):
                            setattr(self,parname+'flag',np.int(lsplit[2]))

                            if (np.int(lsplit[2]) == 1):
                                self.fit_parameters.append(parname)

                        elif (parname != 'START' and parname != 'FINISH' and not is_flag):
                            setattr(self,parname+'err',np.float(lsplit[2]))

                    if (len(lsplit) > 3):
                        if (lsplit[3].find('D') != -1):
                            lsplit[3] = lsplit[3].replace('D','e')
                        setattr(self,parname+'err',np.float(efac*lsplit[3]))

                # store JUMP/EFAC/EQUAD as float, but values have different indeces.
                elif (parname in parameter_list_error):
                    if (parname == 'JUMP'):
                        parname += '_' + lsplit[2]
                        self.fit_parameters.append(parname)
                        setattr(self, parname, np.float(lsplit[3]))
                        if (len(lsplit) > 4):
                            setattr(self, parname + 'flag', efac * np.float(lsplit[4]))
                        if (len(lsplit) > 5):
                            setattr(self, parname + 'err', efac * np.float(lsplit[5]))
                
                    # treat TEMPO2 parameter names somewhat differently.
                    elif (parname is not "JUMP"):
                        setattr(self, "{0}_{1}".format(parname, lsplit[2]) , np.float(lsplit[3]))

                    else:
                        setattr(self, parname, np.float(parvalue))
                        if (len(lsplit) > 2):
                            if (np.int(lsplit[2] == 1)):
                                self.fit_parameters.append(parname)
                            setattr(self,parname+'flag',np.int(lsplit[2]))
                        if (len(lsplit) > 3):
                            setattr(self,parname+'err',np.float(efac*lsplit[3]))

                parorder.append(parname)

            # set array of ordered parameters as object attribute.
            setattr(self, 'parorder', parorder)

    def fix(self):
        """
        Fixes parameters by changing flags to 0.
        """

        for parameter in self.parorder:
            if (hasattr(self, parameter + 'flag')):
                setattr(self, parameter + 'flag', 0)

    def rotate(self, new_epoch, fix_T0=False, rotate_binary_to_new_epoch=False):
        """
        Rotates spin/binary parameters to new PEPOCH, if time-derivatives are present. Default is 
        to rotate solution to the midpoint of the timespan.

        Input
        -----

        new_epoch : float
            New reference epoch for timing solution.

        fix_T0 : bool
            If set to True, leave binary T0 fixed to original value.

        rotate_binary_to_new_epoch : bool
            If set to True, then rotate orbital elements to values exactly at the new PEPOCH.
            This differs from the default operation, where elements are rotated to the rotated 
            T0/TASC value depending on the orbital period. This option is mainly useful for 
            determining the starting elements for numerical integrators. 

        """

        from math import factorial

        old_epoch = self.PEPOCH
        diff_epoch = (new_epoch - old_epoch) * 86400

        # if proper-motion terms are set, rotate sky coordinates.
        if (hasattr(self, 'PMBETA') and hasattr(self, 'PMLAMBDA')):
            pmbeta = self.PMBETA / 1000 / 3600 / 365.25 / 86400
            new_beta = self.BETA + pmbeta * diff_epoch
            setattr(self, 'BETA', new_beta)
            pmlambda = self.PMLAMBDA / 1000 / 3600 / 365.25 / 86400
            new_lambda = self.LAMBDA + pmlambda * diff_epoch
            setattr(self, 'LAMBDA', new_lambda)

        # rotate spin parameters.
        current_spin_location = 0
        for current_spin_1 in parameter_list_spin:
            if hasattr(self, current_spin_1):
                spin_par_1 = getattr(self, current_spin_1)
                idx = 1

                for current_spin_2 in parameter_list_spin[current_spin_location+1:]:
                    fac = factorial(idx)

                    if hasattr(self, current_spin_2):
                        spin_par_2 = getattr(self, current_spin_2)
                        spin_par_1 += spin_par_2 * diff_epoch**idx / fac
                        idx += 1

                    else:
                        break

                setattr(self, 'F' + str(spin_num), spin_par)
            else:
                break

            current_spin_location += 1

        # if in binary system, rotate relevant binary parameters.
        if hasattr(self, 'BINARY'):
            binary_model = getattr(self, 'BINARY')
            new_Tasc = 0.
            old_T0 = self.T0
            new_T0 = 0.

            if (binary_model in model_list_binary_circular):
                pass

            elif (binary_model in model_list_binary_eccentric):
                pb = 0.

                # if BTX, orbital frequencies are used.
                if (binary_model == 'BTX'):
                    pb = 1 / self.FB0 / 86400

                else:
                    pb = self.PB

                # change T0 unless specified otherwise.
                if (not fix_T0):
                    n_orbits = np.int((new_epoch - old_T0) / pb)
                    new_T0 = old_T0 + pb * n_orbits
                    setattr(self, 'T0', new_T0)

                # if desired, set binary time difference to be between 
                # T0 and the new epoch.
                if (rotate_binary_to_new_epoch):
                    new_T0 = new_epoch

            diff_binary = (new_T0 - old_T0) * 86400

            # if derivatives in OM are set, then rotate OM.
            new_OM = self.OM
            idx = 1

            for current_derivative_om in parameter_list_orbit_derivatives["OM"]:
                if hasattr(self, current_derivative_om):
                    current_value = getattr(self, current_derivative_om)
                    fac = factorial(idx)

                    if (current_derivative_om == "OMDOT"):
                        current_value /= (365.25 * 86400) 

                    else:
                        current_value *= (180 / np.pi)

                    new_OM += current_value * diff_binary**idx / fac

            setattr(self, "OM", new_OM)

            # if derivatives in A1 are set, then rotate A1.
            new_A1 = self.A1
            idx = 1

            for current_derivative_a1 in parameter_list_orbit_derivatives["A1"]:
                if hasattr(self, current_derivative_a1):
                    current_value = getattr(self, current_derivative_a1)
                    fac = factorial(idx)

                    if (current_derivative_om == "XDOT"):
                        current_value *= 1e-12

                    new_A1 += current_value * diff_binary**idx / fac

            setattr(self, "A1", new_A1)

            # if derivatves in FB0 are set, then rotate FB0.
            new_FB0 = self.FB0
            idx = 1

            for current_derivative_fb in parameter_list_orbit_derivatives["FB"]:
                if hasattr(self, current_derivative_fb):
                    current_value = getattr(self, current_derivative_fb)
                    fac = factorial(idx)
                    new_FB0 += current_value * diff_binary**idx / fac

            setattr(self, "FB0", new_FB0)

            # if PBDOT is set, then rotate PB.
            if hasattr(self, "PBDOT"):
                pbdot = getattr(self, "PBDOT") * 1e-12
                new_PB = getattr(self, "PB") + (pbdot * diff_binary) / 86400
                setattr(self, "PB", new_PB)
            
    def step(self, uniform_factor=1):
        """
        Randomly steps all parameters with uncertainties in parifle.
        """
        
        n_dim = len(self.fit_parameters)

        for parameter in self.fit_parameters:
            if (hasattr(self,parameter + 'err')):
                value = getattr(self, parameter) 
                err = getattr(self, parameter + 'err') * uniform_factor / n_dim
                
                if (err != 0.):

                    if (parameter == 'RAJ'):
                        ra = Angle(getattr(self, 'RAJ'), unit=u.hour)
                        err = getattr(self, 'RAJerr') / 3600
                        ra_new = ra.deg + err * np.random.uniform(-1., 1.) 
                        ra_new = Angle(ra_new, unit=u.deg)
                        setattr(self, parameter, str(ra_new.to_string(unit=u.hour, sep=':', precision=10)))
    
                    elif (parameter == 'DECJ'):
                        dec = Angle(getattr(self, 'DECJ'), unit=u.deg)
                        err = getattr(self, 'DECJerr') / 3600
                        dec_new = dec.deg + err * np.random.uniform(-1., 1.) 
                        dec_new = Angle(dec_new, unit=u.deg)
                        setattr(self, parameter, str(dec_new.to_string(unit=u.deg, sep=':', precision=10)))

                    elif (parameter == 'SINI'):
                        cosi = np.sqrt(1 - getattr(self, 'SINI')**2)
                        err = getattr(self, 'DECJerr') / 3600
                        cosi +=  err * np.random.uniform(-1., 1.)
                        setattr(self, parameter, np.sqrt(1.-cosi**2))

                    else:
                        setattr(self, parameter, value + err * np.random.uniform(-1., 1.))

                setattr(self,parameter + 'flag', 0)

    def write(self, outfile="out.par"):
        """
        Writes parameter-file Python object to ASCII file.
        """

        from PSRpy.parfile import PrintPar

        PrintPar(self, outfile=outfile)
