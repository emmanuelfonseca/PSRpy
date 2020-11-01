#! /usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys

class DMXPar():
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
