#! /usr/bin/python

import numpy as np
import sys

class ReadTOAs:

    def __init__(self, timfile):

        self.TOAs = []
        self.jump_parameters = []
        self.jump_labels = []
        self.TOA_uncertainties = []

        jump_count = 0
 
        for line in file(timfile):

            if ("MODE" in line):
                pass

            elif ("EFAC" in line):
                pass

            elif ("EQUAD" in line):
                pass

            elif ("FORMAT" in line):
                pass

            elif ("INCLUDE" in line):
                pass

            elif ("JUMP" in line):

                jump_count += 1

                if (jump_count % 2 != 0):
                    if (jump_count < 10):
                        self.jump_parameters.append('JUMP_000' + str(jump_count - (jump_count - 1) / 2))
                    if (jump_count >= 10 and jump_count < 100):
                        self.jump_parameters.append('JUMP_00' + str(jump_count - (jump_count - 1) / 2))

            elif ("INFO" in line):
                pass

            elif ("MODE" in line):
                pass

            else:
                
                elems = line.split()
                self.TOAs.append(np.longdouble(elems[0]))
                self.TOA_uncertainties.append(np.longdouble(elems[1]))
                
                if (jump_count % 2 == 0):
                    self.jump_labels.append('base')

                else:
                    self.jump_labels.append('JUMP_000' + str(jump_count - (jump_count - 1) / 2))

        self.TOAs = np.array(self.TOAs, dtype=np.longdouble)
        self.TOA_uncertainties = np.array(self.TOA_uncertainties, dtype=np.longdouble)
