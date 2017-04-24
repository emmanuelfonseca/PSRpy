#! /usr/bin/python

from re import match
from PSRpy.const import c, G, M_sun, T_sun
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sys

par_strings = ['PSR','PSRJ','RAJ','DECJ','EPHEM','ECL','CLK','UNITS','TIMEEPH',
               'T2CMETHOD','TZRSITE','CORRECT_TROPOSPHERE','PLANET_SHAPIRO',
               'DILATEFREQ','INFO','BINARY','DCOVFILE']
par_ints    = ['NTOA','NITS','NDDM','EPHVER']
par_errors  = ['JUMP', 'T2EFAC', 'T2EQUAD', 'TNECORR', 'ECORR']

pi   = np.pi

class ReadPar():

    def __init__(self, infile, efac=1):
        """
        Reads in parfile parameter values, uncertainties and flags as object attributes.

        Required argument: 
            - 'infile' = input parfile (currently works in TEMPO format).
        Default argument:
            - 'efac' = constant factor to multiply across all uncertainties.
    
        Notes:
            - Parameters stored using same parfile name
              (e.g. par.RAJ stores 'RAJ' value in 'par' object).
            - Parameter uncertainties stored as parameter name plus 'err' 
              extension (e.g. par.RAJerr stores 'RAJ' error in 'par' object).
        """

        # preserve order of parameters in parfile.
        parorder = []

        for line in file(infile):
            lsplit = line.split()

            # if this line is commented out, then skip.
            if (lsplit[0] == 'C' or lsplit[0] == '#'):
                continue

            # else, loop will continue.
            parname, parvalue = lsplit[0], lsplit[1]

            # set the following attributes as strings.
            if (parname in par_strings):
                setattr(self, parname, parvalue)
                # the following is for 'RAJ', 'DECJ' that have flags/errors.
                if (len(lsplit) > 2):
                    setattr(self, parname + 'flag', np.int(lsplit[2]))
                if (len(lsplit) > 3):
                    setattr(self, parname + 'err', efac * lsplit[3])

            # set these as integers.
            elif (parname in par_ints):
                setattr(self, parname, np.int(parvalue))

            # otherwise, if not a JUMP, assume it's a fit parameter 
            # and store value/errors as floats.
            elif (parname not in par_errors):
                # switch 'D' with 'e' for exponents.
                if (parvalue.find('D') != -1):
                    parvalue = parvalue.replace('D','e')
                # DDK model in TEMPO can have SINI = KIN.
                if (parvalue == 'KIN'):
                    setattr(self, parname, parvalue)
                else:
                    setattr(self, parname, np.float(parvalue))
                # set flag and error attributes if present in parfile.
                if (len(lsplit) > 2):
                    setattr(self,parname+'flag',np.int(lsplit[2]))
                if (len(lsplit) > 3):
                    if (lsplit[3].find('D') != -1):
                        lsplit[3] = lsplit[3].replace('D','e')
                    setattr(self,parname+'err',np.float(efac*lsplit[3]))

            # store JUMP/EFAC/EQUAD as float, but values have different indeces.
            elif (parname in par_errors):
                parname += '_'+lsplit[2]
                setattr(self,parname,lsplit[3])
                if (len(lsplit) > 4):
                    setattr(self,parname+'flag',efac*lsplit[4])
                if (len(lsplit) > 5):
                    setattr(self,parname+'err',efac*lsplit[5])

            parorder.append(parname)

        # set array of ordered parameters as object attribute.
        setattr(self,'parorder',parorder)

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
            

    def step(self):
        """
        Randomly steps all parameters with uncertainties in parifle.
        """

        for parameter in self.parorder:
            if (hasattr(self,parameter + 'err')):
                value = getattr(self, parameter) 
                err = getattr(self, parameter + 'err')
                if (parameter == 'SINI'):
                    cosi = np.random.uniform(0.,1.)
                    setattr(self, parameter, np.sqrt(1.-cosi**2))
                else:
                    setattr(self, parameter, value + str(np.random.normal(loc=0., scale=err)))
        setattr(self,parameter+'flag',0)

class DmxPar():
    """
    Extract and manipulate DMX data from parfile object.
    """
    def __init__(self,inobj):
        # check if DMX is set; if so, proceed.
        if (hasattr(inobj,"DMX")):
            nbins   = 0
            dmxgrad = 0
            DMX, DMX1, DMXEP, DMXF1, DMXF2, DMXR1, DMXR2 = [], [], [], [], [], [], []
            DMXerr, DMX1err = [], []
            setattr(self,'DMX',inobj.DMX)
            for ii in range(1000):
                if (ii < 10 and hasattr(inobj,'DMX_000'+str(ii+1))):
                    #print getattr(inobj,'DMX_000'+str(ii+1))
                    DMX.append(getattr(inobj,'DMX_000'+str(ii+1)))
                    DMXerr.append(getattr(inobj,'DMX_000'+str(ii+1)+'err'))
                    if (hasattr(inobj,'DMX1_000'+str(ii+1))):
                       DMX1.append(getattr(inobj,'DMX1_000'+str(ii+1))) 
                       DMX1err.append(getattr(inobj,'DMX1_000'+str(ii+1)+'err'))
                    DMXEP.append(getattr(inobj,'DMXEP_000'+str(ii+1)))
                    DMXF1.append(getattr(inobj,'DMXF1_000'+str(ii+1)))
                    DMXF2.append(getattr(inobj,'DMXF2_000'+str(ii+1)))
                    DMXR1.append(getattr(inobj,'DMXR1_000'+str(ii+1)))
                    DMXR2.append(getattr(inobj,'DMXR2_000'+str(ii+1)))
                elif (ii < 99 and hasattr(inobj,'DMX_00'+str(ii+1))):
                    #print getattr(inobj,'DMX_00'+str(ii+1))
                    DMX.append(getattr(inobj,'DMX_00'+str(ii+1)))
                    DMXerr.append(getattr(inobj,'DMX_00'+str(ii+1)+'err'))
                    if (hasattr(inobj,'DMX1_00'+str(ii+1))):
                       DMX1.append(getattr(inobj,'DMX1_00'+str(ii+1))) 
                       DMX1err.append(getattr(inobj,'DMX1_00'+str(ii+1)+'err'))
                    DMXEP.append(getattr(inobj,'DMXEP_00'+str(ii+1)))
                    DMXF1.append(getattr(inobj,'DMXF1_00'+str(ii+1)))
                    DMXF2.append(getattr(inobj,'DMXF2_00'+str(ii+1)))
                    DMXR1.append(getattr(inobj,'DMXR1_00'+str(ii+1)))
                    DMXR2.append(getattr(inobj,'DMXR2_00'+str(ii+1)))
                elif (ii < 999 and hasattr(inobj,'DMX_0'+str(ii+1))):
                    #print getattr(inobj,'DMX_00'+str(ii+1))
                    DMX.append(getattr(inobj,'DMX_0'+str(ii+1)))
                    DMXerr.append(getattr(inobj,'DMX_0'+str(ii+1)+'err'))
                    if (hasattr(inobj,'DMX1_0'+str(ii+1))):
                       DMX1.append(getattr(inobj,'DMX1_0'+str(ii+1))) 
                       DMX1err.append(getattr(inobj,'DMX1_0'+str(ii+1)+'err'))
                    DMXEP.append(getattr(inobj,'DMXEP_0'+str(ii+1)))
                    DMXF1.append(getattr(inobj,'DMXF1_0'+str(ii+1)))
                    DMXF2.append(getattr(inobj,'DMXF2_0'+str(ii+1)))
                    DMXR1.append(getattr(inobj,'DMXR1_0'+str(ii+1)))
                    DMXR2.append(getattr(inobj,'DMXR2_0'+str(ii+1)))
                else:
                    break
                nbins += 1
        setattr(self,'nbins',nbins)
        setattr(self,'DMX',DMX)
        setattr(self,'DMXerr',DMXerr)
        if (dmxgrad > 0):
            setattr(self,'DMX1',DMX1)
            setattr(self,'DMX1err',DMX1err)
        setattr(self,'DMXEP',DMXEP)
        setattr(self,'DMXR1',DMXR1)
        setattr(self,'DMXR2',DMXR2)
        setattr(self,'DMXF1',DMXF1)
        setattr(self,'DMXF2',DMXF2)
    def plot(self,plotder='n'):
        """
        Plot DMX vs. time with uncertainties (if available).
        """
        plt.errorbar(self.DMXEP,self.DMX,yerr=self.DMXerr,fmt='ro')
        plt.show()


class ConvertPar():
    def __init__(self,inobj,binary=None):
        """
        Converts a subset of data in a supplied parfile object to another 
        desired subset (e.g. converts from equatorial to ecliptic coordinates, or 
        from ELL1 binary model to DD model, if applicable.)
        """
        # if binary option is set...
        if (binary is not None):
            if (binary == inobj.BINARY):
                sys.exit("Error: cannot convert parfile since binary models match!")
            else:
                print "No!"
        #else:
        #    print "No!"
