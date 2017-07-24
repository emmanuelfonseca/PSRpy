#! /usr/bin/python

import numpy as np
import sys

class ReadTOAs:

    def __init__(self, timfile):

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
                print line

            elif ("INFO" in line):
                pass

            elif ("MODE" in line):
                pass

            elif (not line):
                print line

    def stats(self):
        pass
