#! /usr/bin/python

import const
import numpy as np

def doppler(d,b,l,mu):
  """
  Calculate Doppler-bias contribution to orbital decay measurement.
  """
  # header stuff.
  c, pc, d2r = const.c, const.pc, const.d2r
  R0, R0err = np.array([const.R0,const.R0err])*pc*1000.   # convert to m
  v0, v0err = np.array([const.v0,const.v0err])*1000.      # convert to m/s
  z = d*np.sin(b)
  d = d*pc*1000.
  mu = mu/1000./3600.*d2r/86400./365.25
  beta = d/R0*np.cos(b)-np.cos(l)
  # expected bias from Galactic potential, diff. acceleration 
  # and Shklovskii effect.
  galpot = -1.08e-19*(1.25*z/np.sqrt(z**2+0.0324)+0.58*z)*np.sin(b)
  galrel = -np.cos(b)*(v0**2)/R0*(np.cos(l)+beta/(np.sin(l)**2+beta**2))/c
  shklov = (mu**2)*d/c
  # now add it all up.
  #return galpot+galrel+shklov
  return np.array([galpot,galrel,shklov])


def distGR(xpbd,xpbderr,pb,b,l,mu,muerr,nmc=500):
  """
  Determine the distance to pulsar-binary system, assuming GR is correct.
  """
  # header stuff.
  pc, d2r = const.pc, const.d2r
  c = const.c/1000./pc
  dist = np.zeros(nmc)
  R0   = np.random.normal(const.R0,const.R0err,size=(nmc))
  v0   = np.random.normal(const.v0,const.v0err,size=(nmc))/pc
  xpbd = np.random.normal(xpbd,xpbderr,size=(nmc))*1e-12
  mu   = np.random.normal(mu,muerr,size=(nmc))/1000./3600.*d2r/86400./365.25
  pb   = pb*86400.
  # now do the fun stuff, using Newton-Raphson method.
  for i in range(nmc):
    d = 1.
    for j in range(100):
      db = d
      z = d*np.sin(b)
      beta = d/R0[i]*np.cos(b)-np.cos(l)  
      galpot = 1.08e-19*(1.25*z/np.sqrt(z**2+0.0324)+0.58*z)*np.sin(b)
      galrel = np.cos(b)*(v0[i]**2)/R0[i]*(np.cos(l)+beta/(np.sin(l)**2+beta**2))/c
      shklov = (mu[i]**2)*d/c
      # take derivatives of contributions - check these again!
      galpotd = 1.08e-19*(1.25*np.sin(b)/np.sqrt(z**2+0.0324)-\
                1.25*z**2*np.sin(b)/(z**2+0.0324)**(1.5)+0.58)*np.sin(b)
      galreld = (v0[i]**2)*np.cos(b)/c/R0[i]*(1./(np.sin(l)**2+beta**2)-\
                 2.*beta**2*np.cos(b)/R0[i]/(np.sin(l)**2+beta**2)**2)
      shklovd = (mu[i]**2)/c
      f  = xpbd[i]+(galpot+galrel-shklov)*pb
      fp = (galpotd+galreld-shklovd)*pb
      d = d-f/fp
      if (np.fabs(d-db) < 1e-5):
        dist[i] = d
        break
  v0 = v0*pc
  xpbd = xpbd*1e12
  mu = mu*1000.*3600./d2r*86400.*365.25
  return dist, R0, v0, mu, xpbd
