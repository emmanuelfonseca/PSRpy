#! /usr/bin/python 

colors_residuals = ['firebrick', 'royalblue', 'darkmagenta', 'dimgray', 'lightcoral', 'blueviolet', 'brown', 
                    'crimson', 'darkolivegreen', 'lightcoral', 'darkcyan']

matrix_dtype = [('m', '<i4'), ('j', '<i4'), ('paramj', '<a5'), ('gcor', '<f8'), ('sig', '<f8')]

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
