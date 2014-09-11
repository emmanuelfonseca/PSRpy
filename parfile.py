#! /usr/bin/python

from coords import *
import numpy as np
import const

class readpar():
  def __init__(self,infile,efac=1.):
    """
    Reads in parfile parameters and uncertainties as attributes.
    
    Notes:
      (1) Parameters stored using same parfile name
          (e.g. par.RAJ stores 'RAJ' value in 'par' object).
      (2) Parameter uncertainties stored as parameter name plus 'err' 
          extension (e.g. par.RAJerr stores 'RAJ' error in 'par' object).
    """
    for line in file(infile):
      lsplit = line.split()
      # if lines are commented, then skip.
      if (lsplit[0] == 'C' or lsplit[0] == '#'):
        continue
      # set following attributes as strings.
      if (lsplit[0] == 'RAJ' or lsplit[0] == 'DECJ' or lsplit[0] == 'PSRJ' or
      lsplit[0] == 'BINARY' or lsplit[0] == 'PSR' or 
      lsplit[0] == 'EPHEM' or lsplit[0] == 'CLK' or 
      lsplit[0] == 'UNITS' or lsplit[0] == 'TIMEEPH' or 
      lsplit[0] == 'T2CMETHOD' or lsplit[0] == 'CORRECT_TROPOSPHERE' or 
      lsplit[0] == 'PLANET_SHAPIRO' or lsplit[0] == 'DILATEFREQ' or 
      lsplit[0] == 'TZRSITE'):
        setattr(self,lsplit[0],lsplit[1])
        if (len(lsplit) > 2):
          setattr(self,lsplit[0]+'err',np.float(lsplit[3]))
        # if RAJ and DECJ is set, convert to degrees.
        if (lsplit[0] == 'RAJ'):
          setattr(self,lsplit[0]+'deg',ra2deg(lsplit[1]))
        if (lsplit[0] == 'DECJ'):
          setattr(self,lsplit[0]+'deg',dec2deg(lsplit[1]))
      # set these as integers.
      elif (lsplit[0] == 'NTOA' or lsplit[0] == 'NITS'):
        setattr(self,lsplit[0],np.int(lsplit[1]))
      # otherwise, if not a JUMP, store as floats.
      elif (lsplit[0] != 'JUMP'):
        if (lsplit[1].find('D') != -1):
          lsplit[1] = lsplit[1].replace('D','e')
        setattr(self,lsplit[0],np.float(lsplit[1]))
        if (len(lsplit) > 2):
          if (lsplit[3].find('D') != -1):
            lsplit[3] = lsplit[3].replace('D','e')
          setattr(self,lsplit[0]+'err',efac*np.float(lsplit[3]))
      # store JUMP as float, but values have different index.
      elif (lsplit[0] == 'JUMP'):
        setattr(self,lsplit[0]+'_'+lsplit[2],np.float(lsplit[3]))
        if (len(lsplit) > 2):
          setattr(self,lsplit[0]+'_'+lsplit[2]+'err',efac*np.float(lsplit[5]))
          

  
class derivepar():
  def __init__(self,inobj):
    """Derive a bunch of quantities from parfile object."""
    # derive period, derivatives and errors.
    if (hasattr(inobj,'F0')):
      setattr(self,'P0',1./inobj.F0)
    if (hasattr(inobj,'F0err')):
      setattr(self,'P0err',inobj.F0err/inobj.F0**2)
    if (hasattr(inobj,'F1')):
      setattr(self,'P1',-inobj.F1/inobj.F0**2)
    if (hasattr(inobj,'F1err')):
      err = np.sqrt((inobj.F1*inobj.F0err/inobj.F0**3)**2+
            (inobj.F1err/inobj.F0**2)**2)
      setattr(self,'P1err',err)
    if (hasattr(inobj,'F2')):
      setattr(self,'P2',-inobj.F2/inobj.F0**2+2.*inobj.F1**2/inobj.F0**2)
    if (hasattr(inobj,'F2err')):
      setattr(self,'P2err',np.sqrt(((2.*inobj.F2/inobj.F0**3-
              6.*inobj.F1/inobj.F0**4)*inobj.F0err)**2+(4.*inobj.F1*inobj.F1err/
              inobj.F0**3)**2+(inobj.F2err/inobj.F0**2)**2))
    if (hasattr(inobj,'F3')):
      setattr(self,'P3',(6.*inobj.F2*inobj.F1*inobj.F0-inobj.F3*inobj.F0**2-
              6.*inobj.F1**3)/inobj.F0**4)
    if (hasattr(inobj,'F3err')):
      setattr(self,'P3err',
      np.sqrt(((2.*inobj.F3/inobj.F0**3-18.*inobj.F2*inobj.F1/inobj.F0**4+
                24.*inobj.F1**3/inobj.F0**5)*inobj.F0err)**2+
              ((6.*inobj.F2/inobj.F0**3-18.*inobj.F1**2/inobj.F0**4)*
                inobj.F1err)**2+(6.*inobj.F2err*inobj.F1/inobj.F0**3)**2+
               (inobj.F3err/inobj.F0**2)**2))
    # derive proper motion.
    if (hasattr(inobj,'PMRA') and hasattr(inobj,'PMDEC')):
      setattr(self,'mu',np.sqrt(inobj.PMRA**2+inobj.PMDEC**2))
      if (hasattr(inobj,'PMRAerr') and hasattr(inobj,'PMDECerr')):
        setattr(self,'muerr',np.sqrt((inobj.PMRA*inobj.PMRAerr)**2+\
        (inobj.PMDEC*inobj.PMDECerr)**2)/self.mu)
    # derive distance from PX.
    if (hasattr(inobj,'PX')):
      setattr(self,'dist',1./inobj.PX)
      if (hasattr(inobj,'PXerr')):
        setattr(self,'disterr',inobj.PXerr/inobj.PX**2)
    # derive Galactic coordinates of source.
    gal_lat, gal_lon = cel2gal(inobj.RAJdeg,inobj.DECJdeg)
    setattr(self,'gal_b',gal_lat)
    setattr(self,'gal_l',gal_lon)
    # derive binary stuff, if applicable.
    if (hasattr(inobj,'BINARY')):
      # get mass function, in units of Msun.
      pb, a1 = inobj.PB*86400., inobj.A1*const.c
      setattr(self,'massfunc',(4.*np.pi**2/const.G*a1**3/pb**2)/const.Msun)
      # get eccentricity, periastron angle for ELL1 model.
      if (inobj.BINARY == 'ELL1'):
        setattr(self,'E',np.sqrt(inobj.EPS1**2+inobj.EPS2**2))
        setattr(self,'OM',np.arctan2(inobj.EPS2,inobj.EPS1)*180./np.pi)


