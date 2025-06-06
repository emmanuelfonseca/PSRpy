#! /usr/bin/python

import numpy as np
import struct
import sys

def read_resid2(resid2_file, info_file=None):
    """
    Reads the resid2.tmp binary file output by TEMPO.
    
    """

    resid2 = open(resid2_file, "rb")
    resid2_list = []
    info_list = []
    toa_count = 0

    # if an into.txt file is supplied, load in TOA labels.
    if (info_file is not None):
        for label in open(info_file, "r").read().splitlines():
            info_list.append(label)

    # loop over entire file and extract TOA data.
    while True:
        try:
            for bytes_to_read in [4, 72, 4]:
                current_bytes = resid2.read(bytes_to_read)

                if (bytes_to_read == 72):
                    resid2_list.append(struct.unpack('9d', current_bytes))

        except:
            break


    # now define list of dtypes for defining the Numpy array to hold data.
    dtype_list = [
        ('toas', float),
        ('residuals_phase', float),
        ('residuals', float),
        ('orbital_phase', float),
        ('frequency', float),
        ('weight', float),
        ('toa_uncertainties', float),
        ('residuals_prefit', float),
        ('residuals_rednoise', float)
    ]

    # load data into a NumPy record array.
    toa_data = np.array(resid2_list, dtype=dtype_list)

    return toa_data, np.array(info_list, dtype=np.str_)
