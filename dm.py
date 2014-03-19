#! /usr/bin/python

import numpy as np
import sys


def calcstruct(inf,nbinsmax=20,rmbias=False):
  """
  Calculate a structure funtion from DMX measurements.

  Required argument:
    - 'inf' = input DMX file (obtained from running 'dmxparse' on a best-fit parfile).

  Default arguments:
    - 'nbinsmax' = maximum number of time-lag bins. (default is 20)
    - 'rmbias' = remove bias term? (default is no/False)
  """
  # read in data, declare arrays.
  mjd, dmx, err = [], [], []
  count = 0

  for line in file(inf):
    line = line.split()
    mjd.append(np.float(line[0]))
    dmx.append(np.float(line[1]))
    err.append(np.float(line[2]))
    count += 1

  dmx, err, mjd = np.array(dmx), np.array(err), np.array(mjd)

  dmxsq    = []  # squared difference in DMX.
  dmxsqerr = []  # error of squared difference.
  rmserr   = []  # sum of squares of DMX errors.
  lags     = []  # time difference ('lag') between DMX pair.

  # calculate all unique squared difference in DMs.
  for i in range(count-2):
    dmb1, dmb2, dmb3 = dmx[i+1:], mjd[i+1:], err[i+1:]
    for j in range(len(dmb1)):
      dmxsq.append((dmx[i]-dmb1[j])**2)
      dmxsqerr.append(2.*np.sqrt(((dmx[i]-dmb1[j])**2)*(err[i]**2+dmb3[j]**2)))
      rmserr.append((err[i])**2+(dmb3[j])**2)
      lags.append(np.fabs(mjd[i]-dmb2[j]))

  dmxsq    = np.array(dmxsq)
  dmxsqerr = np.array(dmxsqerr)
  rmserr   = np.array(rmserr)
  lags     = np.array(lags)

  # print some info.
  minlag, maxlag = np.min(lags), np.max(lags)
  print "Total number of unique values: {}.".format(len(dmxsq))
  print "Min, max lag: {0:.3f}, {1:.3f} days.".format(minlag,maxlag)

  # get lag bins and structure function, evenly spaced in log_10.
  ulag    = np.unique(lags)
  binedge = 10.**(np.linspace(np.log10(minlag-1.),np.log10(maxlag+1.),num=(nbinsmax+1)))
  bins, struct, structerr = [], [], []

  # compute structure function!
  for k in range(nbinsmax):
    ind = np.where((lags >= binedge[k]) & (lags < binedge[k+1]))
    structbias    = 0.
    structbiaserr = 0.
    # if bin is not empty, proceed.
    if ((dmxsq[ind]).size != 0):
      if (rmbias):
        structbias    = np.mean(rmserr[ind])
        structbiaserr = rmserr[ind]
      struct.append(np.mean(dmxsq[ind])-structbias)
      bins.append(np.mean([binedge[k],binedge[k+1]]))
      structerr.append(np.sqrt(np.sum((dmxsqerr[ind])**2+structbiaserr**2))/np.float(dmxsq[ind].size))

  bins       = np.array(bins) 
  struct     = ((2.*np.pi/1400e6/2.41e-16)**2)*np.array(struct)
  structerr  = ((2.*np.pi/1400e6/2.41e-16)**2)*np.array(structerr)

  return bins, struct, structerr
