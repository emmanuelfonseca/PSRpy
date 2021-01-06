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
        self.backends = []
        self.receivers = []
        
        # now read toas.
        self._read_toas(self.timfile)

        # now derive metadata that is of general use.
        self._get_backends()
        self._get_receivers()

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
                elems = current_flag_line.split()
                flag_idx = elems.index("-{0}".format(flag_name))
                flag_data +=[flag_type(elems[flag_idx + 1])]
           
        else:
            print("WARNING: timfile '{0}' apparently has no TOA flags!".format(self.timfile))

        return flag_data

    def summarize(self):
        """
        Prints a summary of the TOA data loaded into this object.
        """

        # print some basic info first.
        print("Summary and statistics on TOAs in {0}:".format(self.timfile))
        print("  * a total of {0} total TOAs for {1} backends".format(
                len(self.toas),
                len(self.backends)
            )
        )
        print("  * there are data from {0} receivers in the timfile".format(
                len(self.receivers)
            )
        )

        # now move on to per-receiver/backend uncertainty data.
        print("  * uncertainty stats for each receiver:")
        all_receiver_data = self.retrieve_flag_data("f")

        for current_receiver in self.receivers:
            current_uncertainties = []

            for current_idx in range(len(all_receiver_data)):
                if (all_receiver_data[current_idx] == current_receiver):
                    current_uncertainties += [self.toa_uncertainties[current_idx]]

            print("    - {0}: {1:.2f}/{2:.2f}/{3:.2f}".format(
                    current_receiver,
                    np.min(current_uncertainties),
                    np.median(current_uncertainties),
                    np.max(current_uncertainties),
                )
            )

        # now print some info on per-receiver/backend timespan data.
        print("  * timespan data for each receiver:")

        for current_receiver in self.receivers:
            current_toas = []

            for current_idx in range(len(all_receiver_data)):
                if (all_receiver_data[current_idx] == current_receiver):
                    current_toas += [self.toas[current_idx]]

            print("    - {0}: {1:.3f}-{2:.3f}".format(
                    current_receiver,
                    np.min(current_toas),
                    np.max(current_toas),
                )
            )
        
        # now print total number of observations per receiver/backend combination.
        print("  * number of scans for each receiver:")

        for current_receiver in self.receivers:
            current_toas = []

            for current_idx in range(len(all_receiver_data)):
                if (all_receiver_data[current_idx] == current_receiver):
                    current_toas += [self.toas[current_idx]]

            print("    - {0}: {1} scans".format(current_receiver, len(current_toas)))

    def _get_backends(self):
        """
        An internal function that retrieves list of unique backends based on TOA flags 
        if they are present.
        """

        try:
            all_backend_data = self.retrieve_flag_data("be")
            self.backends = list(set(all_backend_data))
            del(all_backend_data)

        except:
            print("INFO: no TOA flag for backend name is present.")

    def _get_receivers(self):
        """
        An internal function that retrieves list of unique receivers based on TOA flags 
        if they are present.
        """

        try:
            all_backend_data = self.retrieve_flag_data("f")
            self.receivers = list(set(all_backend_data))
            del(all_backend_data)

        except:
            pass

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

            elif ("TIME" in line):
                pass

            else:                
                elems = line.split()
                starting_idx = 0

                # prior to extracting data, note if TOA is commented out.
                # it's assumed that all TOA lines have the same structure, 
                # but that commented lines have a 'C' or '#' as the first 
                # character, with a whitespace between it and the next datum.
                if (elems[0] == "C" or elems[0] == "#"):
                    pass 

                else:
                    # now parse TOA data, depending on TOA format.
                    if (self.toa_format == "tempo2"):
                        self.toa_files += [str(elems[0])]
                        self.freqs += [np.float(elems[1])]
                        self.toas += [np.float(elems[2])]
                        self.toa_uncertainties += [np.float(elems[3])]
                        self.observatory_codes += [str(elems[4])]
                        self.toa_flags += [" ".join(elems[5:])]

                    else:
                        print("WARNING: need to add code for toa format '{0}'".format(
                                self.toa_format
                            )
                        )

