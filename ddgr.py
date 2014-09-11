#! /usr/bin/python

import const
import numpy as np


def pbdot(m1,m2,pb,e):
  """
  Calculate orbital decay, as expected from GR.
  """
  # declare constant (i.e.non-mass) terms.
  pi, Tsun = np.pi, const.Tsun
  pb = pb*86400.
  fe = (1.+(73./24.*e**2)+(37./96.*e**4))*(1.-e**2)**(-3.5)
  A  = -192.*pi/5.*(pb/2./pi)**(-5./3.)*fe*Tsun**(5./3.)
  # calculate!
  return A*m1*m2*(m1+m2)**(-1./3.)

#def omdot():
#  """
#  Calculate periastron advance, as expected from GR.
#  """


