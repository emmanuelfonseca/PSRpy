#! /usr/bin/python

import numpy as np

def ra2deg(ra):
  """Convert right ascension (HH:MM:SS) to degrees."""
  val = ra.split(':')
  return (np.float(val[0])+np.float(val[1])/60.+np.float(val[2])/3600.)*15.


def dec2deg(dec):
  """Convert declination (DD:MM:SS) to degrees."""
  val = dec.split(':')
  sign = 1.0
  # if declination < 0, then preserve sign.
  if (np.float(val[0]) < 0.):
    sign = -1.0
  return np.float(val[0])+sign*(np.float(val[1])/60.+np.float(val[2])/3600.)

def cel2gal(ra,dec):
  """Converts celestial coordinates (RAJ, DECJ) into Galactic coordinates."""
  d2r = np.pi/180.
  lat = np.arcsin(np.sin(dec*d2r)*np.cos(62.6*d2r)-
      np.cos(dec*d2r)*np.sin((ra-282.25)*d2r)*np.sin(62.6*d2r))
  lon = 33.*d2r-np.arccos(np.cos(dec*d2r)*np.cos((ra-282.25)*d2r)/np.cos(lat))
  return lat, lon

