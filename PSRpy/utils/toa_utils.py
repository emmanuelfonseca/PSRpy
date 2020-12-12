import numpy as np

def readInfo():
    pass

def readResiduals(in_c_file):
    """
    Reads data from the TEMPO 'c.tmp' file and stores them into a single Python dictionary.
    """

    DataDict = {}
    in_mjd, in_res, in_reserr, in_orbphs = [], [], [], []

    for line in open(in_c_file, "r").readlines():
        
        if ('#' not in line):
            elements = line.split()
            
            in_mjd.append(float(elements[6]))
            in_res.append(float(elements[2]))
            in_reserr.append(float(elements[3]))
            in_orbphs.append(float(elements[5]))
            
    # store as dictionary.
    DataDict['mjd'] = np.array(in_mjd)
    DataDict['residuals'] = np.array(in_res)
    DataDict['residuals_err'] = np.array(in_reserr)
    DataDict['orbital_phase'] = np.array(in_orbphs)

    return DataDict
