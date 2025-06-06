#! /usr/bin/python

from PSRpy.const import c, G, M_sun, T_sun
from astropy.coordinates import SkyCoord
from numpy import pi
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sys


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
        elif (hasattr(inobj,'PMBETA') and hasattr(inobj,'PMLAMBDA')):
            setattr(self,'mu',np.sqrt(inobj.PMBETA**2 + inobj.PMLAMBDA**2))
        elif (hasattr(inobj,'PMELAT') and hasattr(inobj,'PMELONG')):
            setattr(self,'mu',np.sqrt(inobj.PMELAT**2 + inobj.PMELONG**2))
        if (hasattr(inobj,'PMRAerr') and hasattr(inobj,'PMDECerr')):
            setattr(self,'muerr',np.sqrt((inobj.PMRA*inobj.PMRAerr)**2+\
            (inobj.PMDEC*inobj.PMDECerr)**2)/self.mu)
        if (hasattr(inobj,'PMBETAerr') and hasattr(inobj,'PMLAMBDAerr')):
            setattr(self,'muerr',np.sqrt((inobj.PMBETA*inobj.PMBETAerr)**2+\
            (inobj.PMLAMBDA*inobj.PMLAMBDAerr)**2)/self.mu)
        # derive distance from PX, if PX is not equal to zero.
        if (hasattr(inobj,'PX')):
            if (float(inobj.PX) != 0.):
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
            if (inobj.BINARY == 'ELL1' or inobj.BINARY == 'ELL1H'):
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
                if hasattr(inobj, 'ECC'):
                    ecc = inobj.ECC
                else:
                    ecc = inobj.E
            # if orthometric SD parameters are used, convert to m2/cosi.
            if hasattr(inobj,'H3'):
                h3 = inobj.H3
                stig = 0.
                if hasattr(inobj, 'STIG'):
                    stig = inobj.STIG
                elif hasattr(inobj, 'VARSIGMA'):
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
                omdot = float(inobj.OMDOT) * np.pi / 180. / 365.25 / 86400.
                mtot_omdot = (omdot * (1 - float(ecc)**2) / (3. * T_sun**(2./3.)) / (float(pb) / (2. * np.pi))**(-5./3.))**(3./2.)
