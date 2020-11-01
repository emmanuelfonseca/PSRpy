#! /usr/bin/python

from re import match
from PSRpy.const import c, G, M_sun, T_sun
from astropy.coordinates import Angle
from .config_parfile import string_list, error_list, int_list
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sys

pi = np.pi

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

class ReadPar():
    """
    Reads in parfile parameter values, uncertainties and flags as object attributes.

    Inputs
    ------

    infile : str
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


    def __init__(self, infile, efac=1):

        # preserve order of parameters in parfile.
        self.fit_parameters = []
        parorder = []

        for line in open(infile, "r").readlines():
            lsplit = line.split()

            if (len(lsplit) != 0):
                # if this line is commented out, then skip.
                if (lsplit[0] == 'C' or lsplit[0] == '#'):
                    continue

                # else, loop will continue.
                parname, parvalue = lsplit[0], lsplit[1]

                # set the following attributes as strings.
                if (parname in string_list):
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
                elif (parname in int_list):
                    setattr(self, parname, np.int(parvalue))

                # otherwise, if not a JUMP, assume it's a fit parameter 
                # and store value/errors as floats.
                elif (parname not in error_list):
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
                elif (parname in error_list):
                    if (parname == 'JUMP'):
                        parname += '_' + lsplit[2]
                        self.fit_parameters.append(parname)
                        setattr(self, parname, np.float(lsplit[3]))
                        if (len(lsplit) > 4):
                            setattr(self, parname + 'flag', efac * np.float(lsplit[4]))
                        if (len(lsplit) > 5):
                            setattr(self, parname + 'err', efac * np.float(lsplit[5]))
                
                    # treat TEMPO2 parameter names somewhat differently.
                    elif (parname in ["T2EFAC", "T2EQUAD", "ECORR"]):
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

    def rotate(self, new_epoch):
        """
        Rotates spin/binary parameters to new PEPOCH, if time-derivatives are present. Default is 
        to rotate solution to the midpoint of the timespan.

        Input
        -----

        new_epoch : float
            New reference epoch for timing solution.
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
        spinfreq_derivatives = ['F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8',
                                'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16',
                                'F17', 'F18', 'F19', 'F20']
        for spin_num in range(len(spinfreq_derivatives)):
            if hasattr(self, 'F' + str(spin_num)):
                spin_par = getattr(self, 'F' + str(spin_num))
                idx = 1

                for spin_num2 in range(spin_num + 1, len(spinfreq_derivatives)):
                    fac = factorial(idx)
                    if hasattr(self, 'F' + str(spin_num2)):
                        spin_par += (getattr(self, 'F' + str(spin_num2))) * diff_epoch**idx / fac
                        idx += 1
                    else:
                        break

                setattr(self, 'F' + str(spin_num), spin_par)
            else:
                break

        # if in binary system, rotate relevant binary parameters.
        if hasattr(self, 'BINARY'):
            circ_models = ['ELL1', 'ELL1H']
            ecc_models = ['DD', 'DDGR', 'BT', 'BTX']
            binary_model = getattr(self, 'BINARY')
            new_Tasc = 0.
            new_T0 = 0.
            diff_binary = 0.

            if (binary_model in circ_models):
                pass

            elif (binary_model in ecc_models):
                pb = 0.

                # if BTX, orbital frequencies are used.
                if (binary_model == 'BTX'):
                    pb = 1 / self.FB0 / 86400
                else:
                    pb = self.PB

                old_T0 = self.T0
                n_orbits = np.int((new_epoch - old_T0) / pb)
                new_T0 = old_T0 + pb * n_orbits
                diff_binary = new_T0 - old_T0
                setattr(self, 'T0', new_T0)

                # if OMDOT is set, rotate OM.
                if hasattr(self, 'OMDOT'):
                    omdot = self.OMDOT  / 365.25 
                    new_OM = self.OM + omdot * diff_binary 
                    setattr(self, 'OM', new_OM)

            # if XDOT is set, rotate A1.
            if hasattr(self, 'XDOT'):
                xdot = self.XDOT * 1e-12
                new_A1 = self.A1 + xdot * diff_binary * 86400
                setattr(self, 'A1', new_A1)

            # if PBDOT set, rotate PB.
            if hasattr(self, 'PBDOT'):
                pbdot = self.PBDOT * 1e-12
                new_PB = self.PB + pbdot * diff_binary 
                setattr(self, 'PB', new_PB)
            
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
