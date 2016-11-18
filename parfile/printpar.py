#! /usr/bin/python

import numpy as np

__all__ = ["printpar"]

paramslist = ["RAJ","DECJ","PMRA","PMDEC","LAMBDA","BETA","PMLAMBDA","PMBETA",
              "PX","F0","F1","F2","F3","F4","F5","F6","F7","F8","F9","F10",
              "F11","F12","F13","F14","F15","F16","F17","F18","F19","F20",
              "DM","DMX","FD1","FD2","FD3","FD4","FD5","DMX","A1","E","T0","PB",
              "OM","OMDOT","GAMMA","PBDOT","XPBDOT","SINI","MTOT","M2","XDOT",
              "DTHETA","DR","JUMP_1"]

stringlist = ["PSR","PEPOCH","START","FINISH","SOLARN0","EPHEM","ECL","CLK",
              "UNITS","TIMEEPH","T2CMETHOD","CORRECT_TROPOSPHERE","PLANET_SHAPIRO",
              "DILATEFREQ","NTOA","TRES","TZRMJD","TZRFRQ","TZRSITE","NITS","INFO",
              "BINARY"]

errlist    = ["JUMP","T2EFAC","T2EQUAD","ECORR"]

class printpar():
    def __init__(self,inobj):
        """
        Prints a Python parfile object to an ASCII file. 
        """
        # open filehand for writing, loop over attributes in object.
        outfile = open("out.par","w")
        for parameter in inobj.parorder:
            value = str(getattr(inobj,parameter))
            if (parameter in stringlist):
                outfile.write("{0:15}       {1:20}\n".format(parameter,value))
            # isolate noise-model parameters.
            elif ('JUMP' in parameter or 'T2EFAC' in parameter or 'T2EQUAD' in parameter or 'ECORR' in parameter):
                par, und, front = parameter.partition('_')
                fe = '-f'
                if (par == 'JUMP'):
                    fe = '-fe'
                if (hasattr(inobj,parameter+"err")):
                    error, flag = str(getattr(inobj,parameter+"err")), getattr(inobj,parameter+"flag")
                    outfile.write("{0:10} {1:5} {2:20} {3:20} {4:10} {5:10}\n".format(par,fe,front,value,flag,str(error)))
                elif (hasattr(inobj,parameter+"flag")):
                    flag = getattr(inobj,parameter+"flag")
                    outfile.write("{0:10} {1:5} {2:20} {3:10} {4:10}\n".format(par,fe,front,value,flag))
                else:
                    outfile.write("{0:10} {1:5} {2:20} {3:20}\n".format(par,fe,front,value))                  
            else: #if (parameter in paramslist):
                if (value.find('E') != -1):
                    value = value.replace('E','D')
                if (hasattr(inobj,parameter+"err")):
                    error, flag = str(getattr(inobj,parameter+"err")), getattr(inobj,parameter+"flag")
                    if (error.find('E') != -1):
                        error = error.replace('E','D')
                    outfile.write("{0:15}    {1:20}  {2:10}    {3:20}\n".format(parameter,value,flag,str(error)))
                elif (hasattr(inobj,parameter+"flag")):
                    flag = getattr(inobj,parameter+"flag")
                    outfile.write("{0:15}    {1:20}  {2:20}\n".format(parameter,value,flag))
                else:
                    outfile.write("{0:15}    {1:20}\n".format(parameter,value))
        outfile.close()
