#! /usr/bin/python

from numpy.core.records import fromfile as ff
import numpy as np
import sys

# externally define several configuration variables used by the Matrix class.
# makes it easier to change.
matrix_dtype = [
    ('m', '<i4'), 
    ('j', '<i4'), 
    ('paramj', '<a5'), 
    ('gcor', '<f8'), 
    ('sig', '<f8')
]

matrix_parameter_labels = {
    '   f0' : 'F0',
    '   f1' : 'F1',
    '   f2' : 'F2',
    ' f3  ' : 'F3',
    ' f4  ' : 'F4',
    ' f5  ' : 'F5',
    ' f6  ' : 'F6',
    ' f7  ' : 'F7',
    ' f8  ' : 'F8',
    ' f9  ' : 'F9',
    ' f*  ' : 'F10',
    '   RA' : 'RAJ',
    '  Dec' : 'DECJ',
    ' pmra' : 'PMRA',
    ' pmdc' : 'PMDEC',
    '    x' : 'A1',
    '   x2' : 'A1_2',
    ' xdot' : 'XDOT',
    ' XD02' : 'XDOT2',
    '    e' : 'E',
    '   e2' : 'E_2',
    ' edot' : 'EDOT',
    '   T0' : 'T0',
    '  T02' : 'T0_2',
    '   Pb' : 'PB',
    '  Pb2' : 'PB_2',
    ' Pbdt' : 'PBDOT',
    ' Xpbd' : 'XPBDOT',
    ' FB00' : 'FB0',
    ' FB01' : 'FB1',
    ' FB02' : 'FB2',
    '   Om' : 'OM',
    '  Om2' : 'OM_2',
    ' Omdt' : 'OMDOT',
    'gamma' : 'GAMMA',
    '   px' : 'PX',
    '   m2' : 'M2',
    '    M' : 'MTOT',
    '    s' : 'SINI',
    'DX001' : 'DMX_0001',
    'D1001' : 'DMX1_0001',
    'DX002' : 'DMX_0002',
    'D1002' : 'DMX1_0002',
    'DX003' : 'DMX_0003',
    'D1003' : 'DMX1_0003',
    'DX004' : 'DMX_0004',
    'D1004' : 'DMX1_0004',
    'DX005' : 'DMX_0005',
    'D1005' : 'DMX1_0005',
    'DX006' : 'DMX_0006',
    'D1006' : 'DMX1_0006',
    ' O001' : 'JUMP_1',
    ' O002' : 'JUMP_2',
    ' O003' : 'JUMP_3'
}

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
        num = int((num - 29) / 8)

        # declare objects to be returned.
        self.parlist = []
        self.covmat  = np.zeros((num, num), dtype=np.float64)
        freqder_count = 0

        # loop to get data from binary file.
        for ii in range(num):
    
            data = ff(fd, dtype=matrix_dtype, shape=1)
            j = int(data['j'][0])
            m = int(data['m'][0])
            current_label = data['paramj'][0].decode('utf-8')

            if (current_label == ' f*  '):
                self.parlist.append('F1'+str(freqder_count))
                freqder_count += 1
            elif ('DX' not in current_label and 'D1' not in current_label):
                self.parlist.append(matrix_parameter_labels[current_label])

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

        print(self.covmat[index1, index2])
