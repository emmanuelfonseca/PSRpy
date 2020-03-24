#! /usr/bin/python

import numpy as np
import struct
import sys

def read_resid2(resid2_file):
    """
    Reads the resid2.tmp binary file output by TEMPO.
    """

    resid2 = open(resid2_file, "rb")
    resid2_list = []

    # loop over entire file and extract TOA data.
    while True:
        try:
            for bytes_to_read in [4, 72, 4]:
                current_bytes = resid2.read(bytes_to_read)

                if (bytes_to_read == 72):
                    resid2_list.append(struct.unpack('9d', current_bytes))

        except:
            break

    # load data into a NumPy record array.
    toa_data = np.array(resid2_list, dtype=[
        ('toas', np.float),
        ('residuals_phase', np.float),
        ('residuals', np.float),
        ('orbital_phase', np.float),
        ('frequency', np.float),
        ('weight', np.float),
        ('toa_uncertainties', np.float),
        ('residuals_prefit', np.float),
        ('residuals_rednoise', np.float)
    ])

    return toa_data
