#! /usr/bin/python

from PSRpy.const import c, G, M_sun, T_sun
from . import config_parfile as config
from astropy.coordinates import Angle
from . import printpar
from re import match
import PSRpy.utils.math as math
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import sys

class Parfile(object):

    def __init__(self):

        # instantiate all possible parfile parameters.
        template_dict = {"value" : None, "flag"  : None, "error" : None}

        for current_parameter in config.parameter_list_full:
            if current_parameter in config.parameter_list_error:
                setattr(self, current_parameter, {}.copy())
                
            else:
                setattr(self, current_parameter, template_dict.copy())        

        # read input data.
        #self._read_parfile(input_parfile, efac=efac)

        # if present, organize DMX data. 
        self.has_dmx = False 
        #self._get_dmx_data()     

        # next, derive parameters based on parfile data.
        # TODO: define internal-use function for this purpose.


    def get_dmx_data(self):
        """
        A private function that organizes DMX data into lists if present in input parfile.

        Returns
        -------
        has_dmx : bool
            a boolean that is set to True if the input parfile has DMX data.

        n_bins_dmx : int
            the number of DMX values in the parfile.
        
        dmx_val : list of floats
            a list of DMX values. 

        dmx_err : list of floats
            a list of the uncertainties for each DMX value.

        dmx_deriv_val : list of floats
            a list of DMX first derivatives.

        dmx_deriv_err : list of floats
            a list of uncertainties for each DMX first derivative.

        dmx_range_freq : list of two-element lists (floats)
            a list of bounds in frequency for each DMX bin.

        dmx_range_tim : list of two-element lists (floats)
            a list of bounds in time for each DMX bin

        dmx_mean_epoch : list of floats
            a list of mean epochs for each DMX bin.
        """

        try:

            assert(hasattr(self, "DMX"))

            # if assertion passes, proceed with reading DMX data.
            self.has_dmx = True
            self.n_bins_dmx = 0
            self.dmx_val = []
            self.dmx_err = []
            self.dmx_deriv_val = []
            self.dmx_deriv_err = []
            self.dmx_range_freq = []
            self.dmx_range_time = []
            self.dmx_mean_epoch = []

            # now loop over all possible labels and extract existing data.
            for ii in range(config.n_bins_DMX):
                current_label = str(ii).zfill(4)

                if (hasattr(self, "DMX_{0}".format(current_label))):
                    self.dmx_val += [getattr(self, "DMX_{0}".format(current_label))]

                if (hasattr(self, "DMX_{0}err".format(current_label))):
                    self.dmx_err += [getattr(self, "DMX_{0}err".format(current_label))]

                if (hasattr(self, "DMX1_{0}".format(current_label))):
                    self.dmx_deriv_val += [getattr(self, "DMX1_{0}".format(current_label))]

                if (hasattr(self, "DMX1_{0}err".format(current_label))):
                    self.dmx_deriv_val += [getattr(self, "DMX1_{0}err".format(current_label))]

                if (hasattr(self, "DMXEP_{0}".format(current_label))):
                    self.dmx_mean_epoch += [getattr(self, "DMXEP_{0}".format(current_label))]

                # store the time/freq range data as nested lists.
                if (hasattr(self, "DMXF1_{0}".format(current_label)) and
                    hasattr(self, "DMXF2_{0}".format(current_label))):
                    self.dmx_range_freq.append(
                        [
                            getattr(self, "DMXF1_{0}".format(current_label)),
                            getattr(self, "DMXF2_{0}".format(current_label))
                        ]
                    )

                if (hasattr(self, "DMXR1_{0}".format(current_label)) and
                    hasattr(self, "DMXR2_{0}".format(current_label))):
                    self.dmx_range_time.append(
                        [
                            getattr(self, "DMXR1_{0}".format(current_label)),
                            getattr(self, "DMXR2_{0}".format(current_label))
                        ]
                    )

            # now set reamining variables.
            self.n_bins_dmx = len(self.dmx_val)

        except:
            print("WARNING: no DMX data round in parfile!")

    def read(self, input_parfile, efac=1):
        """
        Reads in parfile parameter values, uncertainties and flags as object attributes.

        Inputs
        ------

        input_parfile : str
            pulsar-timing parameter file

        efac : float
            constant factor to multiply across all uncertainties.
    
        Notes
        -----
        Parameters stored as class attibutes using same parameter-file name 
        (e.g. par.RAJ stores the 'RAJ' value in 'par' object). Parameter uncertainties 
        are stored as attributes with parameter name plus an 'err' extension 
        (e.g. par.RAJerr stores the 'RAJ' uncertainty in the 'par' object).
        """

        # loop over each parfile line and extract data.
        for line in open(input_parfile, "r").readlines():
            lsplit = line.split()

            if len(lsplit) != 0:
                # if this line is commented out, then skip.
                if lsplit[0] == 'C' or lsplit[0] == '#':
                    continue

                # else, loop will continue, grab current dictionary.
                parname, parvalue = lsplit[0], lsplit[1]
                current_dict = getattr(self, parname)

                # for non-error terms, set the value in the appropriate type.
                if parname in config.parameter_list_string:
                    current_dict["value"] = parvalue

                elif parname in config.parameter_list_int:
                    current_dict["value"] = int(parvalue)
                
                elif parname not in config.parameter_list_error:
                    # for certain TEMPO values in expontential form, 
                    # switch expnonent character from 'D' to 'e'.
                    if (parvalue.find('D') != -1):
                        parvalue = parvalue.replace('D','e')

                    current_dict["value"] = float(parvalue)

                # before moving on to error terms, suss out flag and uncertainty 
                # info if present.
                if parname not in config.parameter_list_error and len(lsplit) > 2:

                    # the following if/else is to treat ATNF parfiles 
                    # that have no fit flags but do have uncertainties.
                    is_flag = math.represents_an_int(lsplit[2])

                    if is_flag:
                        current_dict["flag"] = int(lsplit[2])

                        if len(lsplit) == 4:
                            # for certain TEMPO values in expontential form, 
                            # switch expnonent character from 'D' to 'e'.
                            if (lsplit[3].find('D') != -1):
                                lsplit[3] = lsplit[3].replace('D','e')

                            current_dict["error"] = float(lsplit[3])

                    else:
                        # for certain TEMPO values in expontential form, 
                        # switch expnonent character from 'D' to 'e'.
                        if (lsplit[2].find('D') != -1):
                            lsplit[2] = lsplit[2].replace('D','e')

                        current_dict["flag"] = 0
                        current_dict["error"] = efac * float(lsplit[2])

                # store JUMP/EFAC/EQUAD as float, but values have different indeces.
                elif parname in config.parameter_list_error:
                    if parname == 'JUMP':
                        jump_dict = {
                            "option" : None, "value" : None, "flag" : None, "error" : None
                        }                        

                        jump_dict["option"] = lsplit[1]
                        jump_dict["value"] = np.float(lsplit[3])

                        # if fit flag and error are present, stash those.
                        if len(lsplit) >= 4:
                            jump_dict["flag"] = int(lsplit[4])

                            if len(lsplit) > 5:
                                jump_dict["error"] = efac * float(lsplit[5])
                
                        current_dict[lsplit[2]] = jump_dict

                    else:
                        error_dict = {"option" : lsplit[1], "value" : float(lsplit[3])}
                        current_dict[lsplit[2]] = error_dict

                # finally, stash loaded dictionary to attribute.
                setattr(self, parname, current_dict)

    def fix(self):
        """
        Fixes parameters by changing flags to 0.
        """

        pass

    def rotate(self, new_epoch, fix_T0=False, rotate_binary_to_new_epoch=False):
        """
        Rotates spin/binary parameters to new PEPOCH, if time-derivatives are present. Default is 
        to rotate solution to the midpoint of the timespan.

        Input
        -----

        new_epoch : float
            New reference epoch for timing solution.

        fix_T0 : bool
            If set to True, leave binary T0 fixed to original value.

        rotate_binary_to_new_epoch : bool
            If set to True, then rotate orbital elements to values exactly at the new PEPOCH.
            This differs from the default operation, where elements are rotated to the rotated 
            T0/TASC value depending on the orbital period. This option is mainly useful for 
            determining the starting elements for numerical integrators. 

        """

        from math import factorial

        old_epoch = self.PEPOCH
        diff_epoch = (new_epoch - old_epoch) * 86400

        # if proper-motion terms are set, rotate sky coordinates.
        if (hasattr(self, 'PMBETA') and hasattr(self, 'PMLAMBDA')):
            pmbeta = self.PMBETA / 1000 / 3600 / 365.25 / 86400
            new_beta = self.BETA + pmbeta * diff_epoch
            setattr(self, 'BETA', new_beta)
            pmlambda = self.PMLAMBDA / 1000 / 3600 / 365.25 / 86400
            new_lambda = self.LAMBDA + pmlambda * diff_epoch
            setattr(self, 'LAMBDA', new_lambda)

        # rotate spin parameters.
        current_spin_location = 0
        for current_spin_1 in config.parameter_list_spin:
            if hasattr(self, current_spin_1):
                spin_par_1 = getattr(self, current_spin_1)
                idx = 1

                for current_spin_2 in config.parameter_list_spin[current_spin_location+1:]:
                    fac = factorial(idx)

                    if hasattr(self, current_spin_2):
                        spin_par_2 = getattr(self, current_spin_2)
                        spin_par_1 += spin_par_2 * diff_epoch**idx / fac
                        idx += 1

                    else:
                        break

                setattr(self, 'F' + str(spin_num), spin_par)
            else:
                break

            current_spin_location += 1

        # if in binary system, rotate relevant binary parameters.
        if hasattr(self, 'BINARY'):
            binary_model = getattr(self, 'BINARY')
            new_Tasc = 0.
            old_T0 = self.T0
            new_T0 = 0.

            if (binary_model in model_list_binary_circular):
                pass

            elif (binary_model in model_list_binary_eccentric):
                pb = 0.

                # if BTX, orbital frequencies are used.
                if (binary_model == 'BTX'):
                    pb = 1 / self.FB0 / 86400

                else:
                    pb = self.PB

                # change T0 unless specified otherwise.
                if (not fix_T0):
                    n_orbits = np.int((new_epoch - old_T0) / pb)
                    new_T0 = old_T0 + pb * n_orbits
                    setattr(self, 'T0', new_T0)

                # if desired, set binary time difference to be between 
                # T0 and the new epoch.
                if (rotate_binary_to_new_epoch):
                    new_T0 = new_epoch

            diff_binary = (new_T0 - old_T0) * 86400

            # if derivatives in OM are set, then rotate OM.
            new_OM = self.OM
            idx = 1

            for current_derivative_om in config.parameter_list_orbit_derivatives["OM"]:
                if hasattr(self, current_derivative_om):
                    current_value = getattr(self, current_derivative_om)
                    fac = factorial(idx)

                    if (current_derivative_om == "OMDOT"):
                        current_value /= (365.25 * 86400) 

                    else:
                        current_value *= (180 / np.pi)

                    new_OM += current_value * diff_binary**idx / fac

            setattr(self, "OM", new_OM)

            # if derivatives in A1 are set, then rotate A1.
            new_A1 = self.A1
            idx = 1

            for current_derivative_a1 in config.parameter_list_orbit_derivatives["A1"]:
                if hasattr(self, current_derivative_a1):
                    current_value = getattr(self, current_derivative_a1)
                    fac = factorial(idx)

                    if (current_derivative_om == "XDOT"):
                        current_value *= 1e-12

                    new_A1 += current_value * diff_binary**idx / fac

            setattr(self, "A1", new_A1)

            # if derivatves in FB0 are set, then rotate FB0.
            new_FB0 = self.FB0
            idx = 1

            for current_derivative_fb in config.parameter_list_orbit_derivatives["FB"]:
                if hasattr(self, current_derivative_fb):
                    current_value = getattr(self, current_derivative_fb)
                    fac = factorial(idx)
                    new_FB0 += current_value * diff_binary**idx / fac

            setattr(self, "FB0", new_FB0)

            # if PBDOT is set, then rotate PB.
            if hasattr(self, "PBDOT"):
                pbdot = getattr(self, "PBDOT") * 1e-12
                new_PB = getattr(self, "PB") + (pbdot * diff_binary) / 86400
                setattr(self, "PB", new_PB)
            
    def step(self, uniform_factor=1):
        """
        Randomly steps all parameters with uncertainties in parifle.
        """
        
        n_dim = len(self.fit_parameters)

        for parameter in self.fit_parameters:
            if (hasattr(self,parameter + 'err')):
                value = getattr(self, parameter) 
                err = getattr(self, parameter + 'err') * uniform_factor / n_dim
                
                if (err != 0.):

                    if (parameter == 'RAJ'):
                        ra = Angle(getattr(self, 'RAJ'), unit=u.hour)
                        err = getattr(self, 'RAJerr') / 3600
                        ra_new = ra.deg + err * np.random.uniform(-1., 1.) 
                        ra_new = Angle(ra_new, unit=u.deg)
                        setattr(
                            self, 
                            parameter, 
                            str(
                                ra_new.to_string(
                                    unit=u.hour, 
                                    sep=':', 
                                    precision=10
                                )
                            )
                        )
    
                    elif (parameter == 'DECJ'):
                        dec = Angle(getattr(self, 'DECJ'), unit=u.deg)
                        err = getattr(self, 'DECJerr') / 3600
                        dec_new = dec.deg + err * np.random.uniform(-1., 1.) 
                        dec_new = Angle(dec_new, unit=u.deg)
                        setattr(
                            self, 
                            parameter, 
                            str(
                                dec_new.to_string(
                                    unit=u.deg, 
                                    sep=':', 
                                    precision=10
                                )
                            )
                        )

                    elif (parameter == 'SINI'):
                        cosi = np.sqrt(1 - getattr(self, 'SINI')**2)
                        err = getattr(self, 'DECJerr') / 3600
                        cosi +=  err * np.random.uniform(-1., 1.)
                        setattr(self, parameter, np.sqrt(1 - cosi**2))

                    else:
                        setattr(self, parameter, value + err * np.random.uniform(-1., 1.))

                setattr(self,parameter + 'flag', 0)

    def write(self, outfile="out.par", sig_fig_error=2):
        """
        Writes parameter-file Python object to ASCII file.
        """

        #printpar.PrintPar(self, outfile=outfile)
        file_lines = []

        # the following is new.
        for current_parameter in config.parameter_list_full:
            current_dict = getattr(self, current_parameter)

            # treat error-adjustment and all other parameters separately.
            if current_parameter not in config.parameter_list_error:
                current_value = current_dict["value"]
                current_line = ""
                
                # if attribute is not empty, proceed.
                if current_value is not None:
                    current_flag = current_dict["flag"]

                    # if flag is present, determine line to print depending on 
                    # whether an uncertainty is present or not.
                    if current_flag is not None:
                        current_error = current_dict["error"]

                        if current_error is not None:
                            current_type = "f"
                            length_value = math.order_of_magnitude(current_dict["error"])
                            length_value += sig_fig_error

                            # if parameters are usually presented in exponent form, 
                            # then adjust here as needed.
                            if current_parameter in config.parameter_list_exponent:
                                current_type = "e"
                                length_value = 7

                            length_name = 25 - length_value 
                            current_line = "{0:<{1}} {2:>.{3}{4}}  {5}  {6}\n".format(
                                current_parameter, length_name, current_value, 
                                length_value, current_type, current_flag, current_error
                            )
                            
                    # otherwise, just print parameter name and value
                    else:
                        current_type = "f"
                        length_value = 10
                        length_value += sig_fig_error

                        # if parameters are usually presented in exponent form, 
                        # then adjust here as needed.
                        if current_parameter in config.parameter_list_exponent:
                            current_type = "e"
                            length_value = 7

                        length_name = 25 - length_value 
                        current_line = ""
                        
                        if (current_parameter in config.parameter_list_string or
                            current_parameter in config.parameter_list_int):
                            current_line = "{0:<{1}} {2:>{3}}\n".format(
                                current_parameter, length_name, current_value, length_value
                            )

                        else:
                            current_line = "{0:<{1}} {2:>.{3}}\n".format(
                                current_parameter, length_name, current_value, 
                                length_value, current_type
                            )

                    file_lines.append(current_line)

                else:
                    pass

            # now treat noise and JUMP terms.
            else:

                # if not empty, proceed.
                if current_dict:
                    
                    # loop over dict entries.
                    for current_key in current_dict.keys():
                        current_line = ""
                        current_option = current_dict[current_key]["option"]
                        current_value = current_dict[current_key]["value"]

                        # if not a JUMP, then print out all entries.
                        if current_parameter != "JUMP":
                            current_line = "{0}   {1}   {2}   {3}\n".format(
                                current_parameter, current_option, current_key, 
                                current_value
                            )

                        # else, treat jump as a fit parameter with option.
                        else:
                            current_flag = current_dict[current_key]["flag"]

                            if current_flag is not None:
                                current_error = current_dict[current_key]["error"]
                                current_line = "{0}   {1}   {2}   {3}   {4}".format(
                                    current_parameter, current_option, current_key, 
                                    current_value, current_flag
                                )

                                if current_error is not None:
                                    current_line = "{0}   {1}\n".format(
                                        current_line, current_error
                                    )

                                else:
                                    current_line = "{0}\n".format(current_line)

                            else:
                                current_line = "{0}   {1}   {2}   {3}\n".format(
                                    current_parameter, current_option, current_key, 
                                    current_value
                                )
                                

                        file_lines.append(current_line)

        # finally, write parfile.
        fout = open(outfile, "w")
        fout.writelines(file_lines)
        fout.close()
