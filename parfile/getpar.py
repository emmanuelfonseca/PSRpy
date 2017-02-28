#! /usr/bin/python

from re import match
from ..const import c, G, M_sun, T_sun
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sys

__all__ = ["ReadPar", "DerivePar", "DmxPar", "ConvertPar"]

par_strings = ['PSR','PSRJ','RAJ','DECJ','EPHEM','ECL','CLK','UNITS','TIMEEPH',
               'T2CMETHOD','TZRSITE','CORRECT_TROPOSPHERE','PLANET_SHAPIRO',
               'DILATEFREQ','INFO','BINARY','DCOVFILE']
par_ints    = ['NTOA','NITS','NDDM','EPHVER']
par_errors  = ['JUMP', 'T2EFAC', 'T2EQUAD', 'TNECORR', 'ECORR']

pi   = np.pi

class ReadPar():
    def __init__(self,infile,efac=1):
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
        if (efac != 1):
            efac = efac
        # preserve order of parameters in parfile.
        parorder = []
        for line in file(infile):
            lsplit = line.split()
            # if this line is commented out, then skip.
            if (lsplit[0] == 'C' or lsplit[0] == '#'):
                continue
            parname, parvalue = lsplit[0], lsplit[1]
            # set the following attributes as strings.
            if (parname in par_strings):
                setattr(self,parname,parvalue)
                # the following is for 'RAJ', 'DECJ'
                if (len(lsplit) > 2):
                    setattr(self,parname+'flag',np.int(lsplit[2]))
                if (len(lsplit) > 3):
                    setattr(self,parname+'err',efac*lsplit[3])
            # set these as integers.
            elif (parname in par_ints):
                setattr(self,parname,np.int(parvalue))
            # otherwise, if not a JUMP, assume it's a fit parameter 
            # and store value/errors as floats.
            elif (parname not in par_errors):
                # switch 'D' with 'e' for exponents.
                if (parvalue.find('D') != -1):
                    parvalue = parvalue.replace('D','e')
                if (parvalue == 'KIN'):
                    setattr(self,parname,parvalue)
                else:
                    setattr(self,parname,np.float(parvalue))
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
        setattr(self,'parorder',parorder)
    def fix(self):
        for parameter in self.parorder:
            if (hasattr(self,parameter+'flag')):
                setattr(self,parameter+'flag',0)
    def step(self):
        for parameter in self.parorder:
            if (hasattr(self,parameter+'err')):
                value, err = getattr(self,parameter), getattr(self,parameter+'err')
                if (parameter == 'SINI'):
                    cosi = np.random.uniform(0.,1.)
                    setattr(self,parameter,np.sqrt(1.-cosi**2))
                else:
                    setattr(self,parameter,value+str(np.random.normal(loc=0.,scale=err)))
        setattr(self,parameter+'flag',0)

class DerivePar():

    def __init__(self,inobj):
        """
        Derive various quantities of potential interest from parfile object.

        Required argument:
            - 'inobj' = parfile object of an input parfile.

        Notes:
            - This sticks to "bare bones", i.e. mostly computes parameters that are usually 
              derived and presented in papers as "derived quantities." 
            - This class is useful when trying to convert one parfile to another with different 
              input parameters (e.g. converting 'ELL1' binary model to 'DD' and generating an 
              equivalent parfile with the latter model).
        """
        # derive RA/DEC in degrees, and Galactic coordinates.
        if (hasattr(inobj,'RAJ') and hasattr(inobj,'DECJ')):
            sc = SkyCoord(inobj.RAJ+' '+inobj.DECJ, unit=(u.hourangle, u.deg))
            setattr(self, 'RAJdeg', sc.ra.deg)
            setattr(self, 'DECJdeg', sc.dec.deg)
            setattr(self, 'gal_b', sc.galactic.b.deg)
            setattr(self, 'gal_l', sc.galactic.l.deg)
        # derive period, period derivatives and their errors.
        if (hasattr(inobj,'F0')):
            setattr(self,'P0',1/inobj.F0)
        if (hasattr(inobj,'F0err')):
            setattr(self,'P0err',inobj.F0err/inobj.F0**2)
        if (hasattr(inobj,'F1')):
            setattr(self,'P1',-inobj.F1/inobj.F0**2)
        if (hasattr(inobj,'F1err')):
            err = np.sqrt((inobj.F1*inobj.F0err/inobj.F0**3)**2+
                  (inobj.F1err/inobj.F0**2)**2)
            setattr(self,'P1err',err)
        if (hasattr(inobj,'F2')):
            setattr(self,'P2',-inobj.F2/inobj.F0**2+2*inobj.F1**2/inobj.F0**2)
        if (hasattr(inobj,'F2err')):
            setattr(self,'P2err',np.sqrt(((2*inobj.F2/inobj.F0**3-
                    6*inobj.F1/inobj.F0**4)*inobj.F0err)**2+(4*inobj.F1*inobj.F1err/
                    inobj.F0**3)**2+(inobj.F2err/inobj.F0**2)**2))
        if (hasattr(inobj,'F3')):
            setattr(self,'P3',(6*inobj.F2*inobj.F1*inobj.F0-inobj.F3*inobj.F0**2-
                    6*inobj.F1**3)/inobj.F0**4)
        if (hasattr(inobj,'F3err')):
            setattr(self,'P3err',
            np.sqrt(((2*inobj.F3/inobj.F0**3-18*inobj.F2*inobj.F1/inobj.F0**4+
                      24*inobj.F1**3/inobj.F0**5)*inobj.F0err)**2+
                      ((6*inobj.F2/inobj.F0**3-18*inobj.F1**2/inobj.F0**4)*
                      inobj.F1err)**2+(6*inobj.F2err*inobj.F1/inobj.F0**3)**2+
                      (inobj.F3err/inobj.F0**2)**2))
        # derive proper motion.
        if (hasattr(inobj,'PMRA') and hasattr(inobj,'PMDEC')): 
            setattr(self,'mu',np.sqrt(inobj.PMRA**2 + inobj.PMDEC**2))
        if (hasattr(inobj,'PMBETA') and hasattr(inobj,'PMLAMBDA')):
            setattr(self,'mu',np.sqrt(inobj.PMBETA**2 + inobj.PMLAMBDA**2))
        if (hasattr(inobj,'PMRAerr') and hasattr(inobj,'PMDECerr')):
            setattr(self,'muerr',np.sqrt((inobj.PMRA*inobj.PMRAerr)**2+\
            (inobj.PMDEC*inobj.PMDECerr)**2)/self.mu)
        if (hasattr(inobj,'PMBETAerr') and hasattr(inobj,'PMLAMBDAerr')):
            setattr(self,'muerr',np.sqrt((inobj.PMBETA*inobj.PMBETAerr)**2+\
            (inobj.PMLAMBDA*inobj.PMLAMBDAerr)**2)/self.mu)
        # derive distance from PX, if PX is not equal to zero.
        if (hasattr(inobj,'PX')):
            if (np.float(inobj.PX) != 0.):
                setattr(self,'dist',1/inobj.PX)
            if (hasattr(inobj,'PXerr')):
                setattr(self,'disterr',inobj.PXerr/inobj.PX**2)
        # derive binary stuff, if applicable.
        if (hasattr(inobj,'BINARY')):
            pb, ecc = 0, 0
            if (inobj.BINARY == 'BTX'):
                setattr(self,'PB',1/inobj.FB)
                pb = 1/inobj.FB*86400.
            else:
                pb = inobj.PB*86400
            # if using ELL1 mode, get eccentricity, argument/epoch of periastron.
            if (inobj.BINARY == 'ELL1'):
                eps1, eps2 = inobj.EPS1, inobj.EPS2
                ecc = np.sqrt(eps1**2 + eps2**2)
                setattr(self, 'E', ecc)
                om = np.arctan2(float(inobj.EPS1),float(inobj.EPS2))
                setattr(self,'OM',om*180/pi)
                setattr(self,'T0',inobj.TASC+pb/86400/2/pi*om)
                if (hasattr(inobj, 'EPS1err') and hasattr(inobj, 'EPS2err')):
                    eps1_err, eps2_err = inobj.EPS1err, inobj.EPS2err
                    de = np.sqrt((eps1 * eps1_err)**2 + (eps2 * eps2_err)**2) / ecc
                    setattr(self, 'Eerr', de)
                    dom = np.sqrt((eps2 * eps1_err)**2 + (eps1 * eps2_err)**2) / ecc**2
                    setattr(self, 'OMerr', dom * 180 / pi)
                    Tasc_err = inobj.TASCerr
                    dT0 = np.sqrt(Tasc_err**2 + (inobj.PB / 2*np.pi * dom)**2)
                    setattr(self, 'T0err', dT0)
            else:
                ecc = inobj.E
            # if orthometric SD parameters are used, convert to m2/cosi.
            if hasattr(inobj,'H3'):
                h3 = inobj.H3
                stig = 0.
                if hasattr(inobj, 'STIG'):
                    stig = inobj.STIG
                elif hasattr(inob, 'VARSIGMA'):
                    stig = inobj.VARSIGMA
                m2 = h3/stig**3/T_sun
                cosi = (1-stig**2)/(1+stig**2)
                sini = np.sqrt(1-cosi**2)
                setattr(self,'M2',m2)
                setattr(self,'SINI',sini)
            # get mass function, in units of M_sun.
            a1 = inobj.A1*c 
            mfunc = (4*pi**2/G*a1**3/pb**2)/M_sun     
            setattr(self,'massfunc',mfunc)
            # compute total mass from mass function, if SHapiro delay is measured.
            if (hasattr(inobj, 'M2') and hasattr(inobj, 'SINI')):
                m2, sini = inobj.M2, inobj.SINI
                mtot = np.sqrt((m2 * sini)**3 / mfunc)
                setattr(self, 'MTOT', mtot)
            # compute total mass from OMDOT, if it's there.
            # this assumes that the OMDOT is due to GR.
            if (hasattr(inobj, 'OMDOT') and getattr(inobj, 'OMDOT') > 0):
                omdot = np.float(inobj.OMDOT) * np.pi / 180. / 365.25 / 86400.
                mtot_omdot = (omdot * (1 - np.float(ecc)**2) / (3. * T_sun**(2./3.)) / (np.float(pb) / (2. * np.pi))**(-5./3.))**(3./2.)
                
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
