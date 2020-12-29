#! /usr/bin/python

import numpy as np
import sys

class Timfile(object):

    def __init__(self, timfile, toa_format):

        self.timfile = timfile
        self.toa_format = toa_format
        self.toa_files = []
        self.freqs = []
        self.toas = []
        self.jump_parameters = []
        self.jump_labels = []
        self.observatory_codes = []
        self.toa_uncertainties = []
        self.toa_flags = []
        self.is_toa_removed = []
        
        # now read toas.
        self._read_toas(self.timfile)

    def retrieve_flag_data(self, flag_name, flag_type=str):
        """
        Loops over existing TOA-flag data and extracts the value if flag is present.

        Parameters
        ----------
        flag_name : str
            name of the TOA flag for which to extract data.

        flag_tyoe : type, optional
            a Python type of the desired data (default: str)

        Returns
        -------
        flag_data : list of type 'flag_type'
            a list containing the flag data.
        """

        flag_data = []

        if (len(self.toa_flags) > 0):

            for current_flag_line in self.toa_flags:
                print(current_flag_line)
                elems = current_flag_line.split()
                flag_idx = elems.index("-{0}".format(flag_name))
                flag_data +=[flag_type(elems[flag_idx + 1])]
           
        else:
            print("WARNING: timfile '{0}' apparently has no TOA flags!".format(self.timfile))

        return flag_data

    def _read_toas(self, timfile):

        jump_count = 0
 
        for line in open(timfile, "r"):

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
                pass

            elif ("INFO" in line):
                pass

            elif ("MODE" in line):
                pass

            else:                
                elems = line.split()
                starting_idx = 0

                # prior to extracting data, note if TOA is commented out.
                # it's assumed that all TOA lines have the same structure, 
                # but that commented lines have a 'C' or '#' as the first 
                # character, with a whitespace between it and the next data.
                if (elems[0] == "C" or elems[0] == "#"):
                    self.is_toa_removed += [True]
                    starting_idx = 1

                else:
                    self.is_toa_removed += [False]
                
                # now parse TOA data, depending on TOA format.
                if (self.toa_format == "tempo2"):
                    self.toa_files += [str(elems[0 + starting_idx])]
                    self.freqs += [np.float(elems[1 + starting_idx])]
                    self.toas += [np.float(elems[2 + starting_idx])]
                    self.toa_uncertainties += [np.float(elems[3 + starting_idx])]
                    self.observatory_codes += [str(elems[4 + starting_idx])]
                    self.toa_flags += [" ".join(elems[5 + starting_idx:])]

                else:
                    print("WARNING: need to add code for toa format '{0}'".format(self.toa_format))
