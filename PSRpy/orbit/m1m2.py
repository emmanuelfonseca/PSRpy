#! /usr/bin/python

from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from ..const import T_sun, d2r
from ..parfile import DerivePar
from .pkcorr import doppler
import matplotlib.pyplot as plt
import numpy as np

rcParams["mathtext.fontset"] = "cm"

def om1dot_m1m2(omd,omderr,a1,e,pb,om,m1,npts):
    """
    Calculate upper/lower bounds of OMDOT curve in the m1-m2 plane.
    """

    m2omh = ((omd+omderr)*d2r*(1.-e**2)*(pb/2./np.pi)**(5./3.)/3./\
              T_sun.value**(2./3.)/86400./365.25)**(3./2.) - m1
    m2oml = ((omd-omderr)*d2r*(1.-e**2)*(pb/2./np.pi)**(5./3.)/3./\
              T_sun.value**(2./3.)/86400./365.25)**(3./2.) - m1
    return m2omh, m2oml


def pbdot_m1m2(pbdot,pbdoterr,pb,e,m1,npts):
    """
    Calculate the upper/lower bounds of PBDOT curve in the m1-m2 plane.
    """

    m2pbdh = np.zeros(npts)
    m2pbdl = np.zeros(npts)
    fe = 1.+73./24.*e**2+37./96.*e**4
    A  = -192.*np.pi/5.*(pb/2./np.pi)**(-5./3.)*fe*(1.-e**2)**(-7./2.)*\
         T_sun.value**(5./3.)

    for i in range(npts):
        m2 = 1.
        # use Newton-Raphson method to get upper-bound curve.
        if (m1[i] == 0.):
            m1[i] = 0.001
        for j in range(100):
            m2b = m2
            f   = A*m1[i]*m2*(m1[i]+m2)**(-1./3.)
            fp  = A*m1[i]*m2*(m1[i]+m2)**(-1./3.)*(1./m2-1./3./(m1[i]+m2))
            m2  = m2-((pbdot+pbdoterr)*1e-12-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2pbdh[i] = m2
                break

    for k in range(npts):
        m2 = 1.
        # use Newton-Raphson method to get upper-bound curve.
        if (m1[k] == 0.):
            m1[k] = 0.001
        for l in range(100):
            m2b = m2
            f   = A*m1[k]*m2*(m1[k]+m2)**(-1./3.)
            fp  = A*m1[k]*m2*(m1[k]+m2)**(-1./3.)*(1./m2-1./3./(m1[k]+m2))
            m2  = m2-((pbdot-pbdoterr)*1e-12-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2pbdl[k] = m2
                break

    return m2pbdh, m2pbdl


def gamma_m1m2(gam,gamerr,e,pb,m1,npts):
    """
    Calculate upper/lower bounds of GAMMA curve in the m1-m2 plane.
    """

    m2gamh = np.zeros(npts)
    m2gaml = np.zeros(npts)

    for i in range(npts):
        m2 = 1.
        # use Newton-Raphson method to get upper-bound curve.
        for j in range(100):
            m2b = m2
            f   = e*(pb/2./np.pi)**(1./3.)*T_sun.value**(2./3.)*(m1[i]+m2)**(-4./3.)*\
                  m2*(m1[i]+2.*m2)
            fp  = e*(pb/2./np.pi)**(1./3.)*T_sun.value**(2./3.)*(m1[i]+m2)**(-4./3.)*\
                  (-4./3./(m1[i]+m2)*m2*(m1[i]+m2)+m1[i]+4.*m2) 
            m2  = m2-(gam+gamerr-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2gamh[i] = m2
                break

    for k in range(npts):
        m2 = 1.
        # use Newton-Raphson method to get lower-bound curve.
        for l in range(100):
            m2b = m2
            f   = e*(pb/2./np.pi)**(1./3.)*T_sun.value**(2./3.)*(m1[k]+m2)**(-4./3.)*\
                  m2*(m1[k]+2.*m2)
            fp  = e*(pb/2./np.pi)**(1./3.)*T_sun.value**(2./3.)*(m1[k]+m2)**(-4./3.)*\
                  (-4./3./(m1[k]+m2)*m2*(m1[k]+2.*m2)+m1[k]+4.*m2) 
            m2  = m2-(gam-gamerr-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2gaml[k] = m2
                break

    return m2gamh, m2gaml


def r_m1m2(m2,m2err,npts):
    """
    Calculate upper/lower bounds of Shapiro-r curve in the m1-m2 plane.
    """

    m2rh = np.zeros(npts)+(m2+m2err)
    m2rl = np.zeros(npts)+(m2-m2err)

    return m2rh, m2rl


def s_m1m2(s,serr,x,pb,m1,npts):
    """
    Calculate upper/lower bounds of Shapiro-s curve in the m1-m2 plane.
    """

    m2sh = np.zeros(npts)
    m2sl = np.zeros(npts)

    for i in range(npts):
        m2 = 0.2
        # use Newton-Raphson method to get upper-bound curve.
        for j in range(100):
            m2b = m2
            f   = x*(pb/2./np.pi)**(-2./3.)*T_sun.value**(-1./3.)*(m1[i]+m2)**(2./3.)/m2
            fp  = x*(pb/2./np.pi)**(-2./3.)*T_sun.value**(-1./3.)*(m1[i]+m2)**(2./3.)/m2*\
                  (2./3./(m1[i]+m2)-1./m2)
            m2  = m2-(s+serr-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2sh[i] = m2
                break

    for k in range(npts):
        m2 = 0.2
        # use Newton-Raphson method to get upper-bound curve.
        for l in range(100):
            m2b = m2
            f   = x*(pb/2./np.pi)**(-2./3.)*T_sun.value**(-1./3.)*(m1[k]+m2)**(2./3.)/m2
            fp  = x*(pb/2./np.pi)**(-2./3.)*T_sun.value**(-1./3.)*(m1[k]+m2)**(2./3.)/m2*\
                  (2./3./(m1[k]+m2)-1./m2)
            m2  = m2-(s-serr-f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2sl[k] = m2
                break

    return m2sh, m2sl

def om1spin_m1m2(om1s,om1serr,pb,ecc,m1,npts):
    """
    Calculate upper/lower bounds of precession-rate curve in the m1-m2 plane.
    """

    om1serr_lo, om1serr_up = om1serr[0], om1serr[1]
    m2sh = np.zeros(npts)
    m2sl = np.zeros(npts)
    A = 0.5 * (T_sun.value)**(2./3.) * (pb / 2 / np.pi)**(-5./3.) / \
        (1 - ecc**2) * 86400 * 365.25

    for i in range(npts):
        m2 = 0.2
        # use Newton-Raphson method to get upper-bound curve.
        for j in range(100):
            m2b = m2
            f  = A * (m2 * (4*m1[i] + 3*m2) / (m1[i] + m2)**(4./3.)) / d2r
            fp = A * ((m1[i] + m2) * (4*m1[i] + 6*m2) - (4./3.) * m2 * (4*m1[i] + 3*m2)) / \
                 (m1[i] + m2)**(7./3.) / d2r
            m2 = m2-(om1s + om1serr_up - f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2sh[i] = m2
                break

    for k in range(npts):
        m2 = 0.2
        for l in range(100):
            m2b = m2
            f  = A * (m2 * (4*m1[k] + 3*m2) / (m1[k] + m2)**(4./3.)) / d2r
            fp = A * ((m1[k] + m2) * (4*m1[k] + 6*m2) - (4./3.) * m2 * (4*m1[k] + 3*m2)) / \
                 (m1[k] + m2)**(7./3.) / d2r
            m2 = m2 - (om1s - om1serr_lo - f)/(-fp)
            if (np.fabs(m2-m2b) < 1e-7):
                m2sl[k] = m2
                break

    return m2sh, m2sl

# define m1m2 class that uses the above functions when they're set in parfile.

class M1M2(object):

    def __init__(self, inobj, npts=200, font_name="serif", om1s=[False,False,False], pkcorr='n'):
        """
        Calculate upper/lower bounds of post-Keplerian (PK)parameters and store 
        all arrays in a single object. Currently supported PK parameters are: 
        PBDOT, OMDOT, GAMMA, r, s.

        Required argument:
            - 'inobj'  = input parfile object, which can be generated by using 
                         the "readpar" class in the 'psrpar.py' module.
        Default arguments:
            - 'npts'   = number of array elements for upper/lower bounds.
            - 'pkcorr' = correct for Doppler bias? (y = yes, n = no).
            - 'om1s'   = list of values for a geodetic-precession measurement.
                         (median value, lower uncertainy, upper uncertainy.)
        """

        self.font = FontProperties()
        self.font.set_name(font_name)

        a1 = inobj.A1["value"] 
        e = inobj.E["value"] 
        pb = inobj.PB["value"] * 86400
        om = inobj.OM["value"]
        m1 = 3. * np.arange(npts) / (npts-1.)

        setattr(self,'m1',m1)
        m2omh = m2oml = np.zeros(npts)

        # if OMDOT and its error are set, calculate m1m2 arrays.
        if inobj.OMDOT["error"] is not None:
            omdot = inobj.OMDOT["value"] 
            omdoterr = inobj.OMDOT["error"]
            m2omh, m2oml = om1dot_m1m2(omdot,omdoterr,a1,e,pb,om,m1,npts)
            setattr(self,'OMDOT_U',m2omh)
            setattr(self,'OMDOT_L',m2oml)

        # if GAMMA and its error are set, calculate m1m2 arrays.
        if inobj.GAMMA["error"] is not None:
            gamma = inobj.GAMMA["value"]
            gammaerr = inobj.GAMMA["error"]
            m2gamh, m2gaml  = gamma_m1m2(gamma,gammaerr,e,pb,m1,npts)
            setattr(self,'GAMMA_U',m2gamh)
            setattr(self,'GAMMA_L',m2gaml)

        # if Shapiro-r and its error are set, calculate m1m2 arrays.
        if inobj.M2["error"] is not None :
            r = inobj.M2["value"] 
            rerr = inobj.M2["error"]
            m2rh, m2rl  = r_m1m2(r,rerr,npts)
            setattr(self,'r_U',m2rh)
            setattr(self,'r_L',m2rl)

        # if Shapiro-s and its error are set, calculate m1m2 arrays.
        if inobj.SINI["error"] is not None:
            s = inobj.SINI["value"] 
            serr = inobj.SINI["error"]
            m2sh, m2sl  = s_m1m2(s,serr,a1,pb,m1,npts)
            setattr(self,'s_U',m2sh)
            setattr(self,'s_L',m2sl)

        # if PBDOT and its error are set, calculate m1m2 arrays.
        if inobj.PBDOT["error"] is not None:
            pbdot = inobj.PBDOT["value"]
            pbdoterr = inobj.PBDOT["error"]
            #if (pkcorr == 'y'):
            #    der = DerivePar(inobj)
            #    b, l, mu, muerr = der.gal_b, der.gal_l, der.mu, der.muerr
            #    corr, corr_err = doppler(0.7,0.3,b,l,mu,muerr)
            #    corr     *= inobj.PB*86400.*1e12
            #    corr_err *= inobj.PB*86400.*1e12
            #    pbdot -= sum(corr)
            #    pbdoterr = np.sqrt(pbdoterr**2+corr_err**2)
            m2pbdh, m2pbdl  = pbdot_m1m2(pbdot,pbdoterr,pb,e,m1,npts)
            setattr(self,'PBDOT_U',m2pbdh)
            setattr(self,'PBDOT_L',m2pbdl)

        # if OM1SPIN and the uncertainties are set, calculate m1m2 arrays.
        if (any(om1s)):
            om1smed = om1s[0]
            om1serr = om1s[1:]
            print(om1serr)
            m2om1sh, m2om1sl = om1spin_m1m2(om1smed, om1serr, pb, e, m1, npts)
            setattr(self,'OM1SPIN_U',m2om1sh)
            setattr(self,'OM1SPIN_L',m2om1sl)

    def plot(self, pbdot_extra: float = None):
        """
        Plot availabe m1-m2 data using the matplotlib package.
        """

        if (hasattr(self,'OMDOT_U') and hasattr(self,'OMDOT_L')):
            #plt.plot(self.m1,self.OMDOT_U,'k-')
            #plt.plot(self.m1,self.OMDOT_L,'k-')
            plt.fill_between(self.m1, self.OMDOT_U, self.OMDOT_L, color='k', alpha=0.8)
            plt.text(2.0, 0.4, r'$\dot{\omega}$', fontproperties=self.font, fontsize=15)

        if (hasattr(self,'GAMMA_U') and hasattr(self,'GAMMA_L')):
            #plt.plot(self.m1,self.GAMMA_U,'g-')
            #plt.plot(self.m1,self.GAMMA_L,'g-')
            plt.text(0.2, 1.0, r'$\gamma$', fontproperties=self.font, fontsize=15)
            plt.fill_between(self.m1, self.GAMMA_U, self.GAMMA_L, color='g', alpha=0.8)

        if (hasattr(self,'r_U') and hasattr(self,'r_L')):
            #plt.plot(self.m1,self.r_U,'r-')
            #plt.plot(self.m1,self.r_L,'r-')
            plt.fill_between(self.m1, self.r_U, self.r_L, color='r', alpha=0.5)
            plt.text(2.7, 1.1, r'$r$', fontproperties=self.font, fontsize=15)

        if (hasattr(self,'s_U') and hasattr(self,'s_L')):
            #plt.plot(self.m1,self.s_U,'m-')
            #plt.plot(self.m1,self.s_L,'m-')
            plt.fill_between(self.m1, self.s_U, self.s_L, color='m', alpha=0.5)
            plt.text(0.5, 0.7, r'$s$', fontproperties=self.font, fontsize=15)

        if (hasattr(self,'PBDOT_U') and hasattr(self,'PBDOT_L')):
            #plt.plot(self.m1,self.PBDOT_U,'b-')
            #plt.plot(self.m1,self.PBDOT_L,'b-')
            plt.fill_between(self.m1, self.PBDOT_U, self.PBDOT_L, color='blue', alpha=0.5)
            plt.text(0.9, 2.5, r'$\dot{P}_{\rm b}$', fontproperties=self.font, fontsize=15)

        if pbdot_extra is not None:
            pbdot_extra_lo, pbdot_extra_hi = pbdot_extra
            plt.plot(self.m1, pbdot_extra_lo, 'k--')
            plt.plot(self.m1, pbdot_extra_hi, 'k--')

        if (hasattr(self,'OM1SPIN_U') and hasattr(self,'OM1SPIN_L')):
            #plt.plot(self.m1,self.OM1SPIN_U,'y-')
            #plt.plot(self.m1,self.OM1SPIN_L,'y-')
            plt.fill_between(self.m1, self.OM1SPIN_U, self.OM1SPIN_L, color='yellow', alpha=0.7)
            plt.text(2.5, 2.3, r'$\Omega_1^{\rm spin}$', fontproperties=self.font, fontsize=15)

        plt.xlim(0.,3.)
        plt.ylim(0.,3.)
        #plt.axes().set_aspect('equal')
        plt.xlabel(r'Pulsar Mass (${\rm M}_{\odot}$)', fontproperties=self.font, fontsize=15)
        plt.ylabel(r'Companion Mass (${\rm M}_{\odot}$)', fontproperties=self.font, fontsize=15)
        plt.savefig('m1m2.png', dpi=500, fmt='png')
        plt.show()

