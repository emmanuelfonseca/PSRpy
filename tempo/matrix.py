#! /usr/bin/python

from config_tempo import matrix_dtype, matrix_parameter_labels
from numpy.core.records import fromfile as ff
import numpy as np
import sys

class Matrix: 
    """
    A Python class that reads and extracts information from a covariance matrix 
    estimated by TEMPO.
    """

    def __init__(self, matrix_file='matrix.tmp'):
    
         # open binary file and get number of params.
        fd = open(matrix_file, 'rb')

        num = ff(fd, formats='i4', shape=1)
        num = int(num[0][0])
        num = (num - 29) / 8

        # declare objects to be returned.
        self.parlist = []
        self.covmat  = np.zeros([num, num], dtype=np.float64)
        freqder_count = 0

        # loop to get data from binary file.
        for ii in range(num):
    
            data = ff(fd, dtype=matrix_dtype, shape=1)
            j = int(data['j'][0])
            m = int(data['m'][0])

            if (data['paramj'][0] == ' f*  '):
                self.parlist.append('F1'+str(freqder_count))
                freqder_count += 1
            else:
                self.parlist.append(matrix_parameter_labels[data['paramj'][0]])

            a = ff(fd, formats='f8', shape=m, byteorder='<')
    
            for jj in range(num):
    
                self.covmat[ii, jj] = a[jj][0]

            fd.seek(8, 1)

    def correlation(self, parameter1, parameter2):
        """
        Returns the covariance between two parameters, using resultant TEMPO matrix file.
        """

        index1 = self.parlist.index(parameter1)
        index2 = self.parlist.index(parameter2)

        print self.covmat[index1, index2]
