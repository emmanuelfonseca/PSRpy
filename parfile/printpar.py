#! /usr/bin/python

from config_parfile import parameter_list, string_list, int_list, error_list, spin_freq_derivatives
import numpy as np
import sys

class PrintPar():
    """
    A class that prints a Parfile object to an ASCII file using the tempo format.
    
    Input
    -----

    inobj : object
        Python object that has parameter-file content stored as attributes.

    outfule : str
        name of output file.
    """

    def __init__(self, inobj, outfile='out.par'):
        """
        Prints a Python parfile object to an ASCII file. 
        """

        outfile = open(outfile, "w")

        for parameter in inobj.parorder:

            value = str(getattr(inobj,parameter))

            # first check if parameter has fixed/string value.
            if (parameter in string_list or parameter in int_list):
                outfile.write("{0:15}       {1:20}\n".format(parameter, value))

            # next, isolate noise-model parameters.
            elif ('JUMP' in parameter or 'T2EFAC' in parameter or 'T2EQUAD' in parameter or 'ECORR' in parameter):
                par, und, front = parameter.partition('_')
                fe = '-f'

                if (par == 'JUMP'):
                    fe = '-fe'

                if (hasattr(inobj,parameter+"err")):
                    error, flag = str(getattr(inobj,parameter+"err")), getattr(inobj,parameter+"flag")
                    outfile.write("{0:10} {1:5} {2:20} {3:20} {4:10} {5:10}\n".format(par,fe,front,value,flag,str(error)))

                elif (hasattr(inobj,parameter+"flag")):
                    flag = getattr(inobj,parameter+"flag")
                    outfile.write("{0:10} {1:5} {2:20} {3:10} {4:10}\n".format(par,fe,front,value,flag))

                else:
                    outfile.write("{0:10} {1:5} {2:20} {3:20}\n".format(par,fe,front,value))                  

            # finally, loop over fit parameters.
            elif (parameter in parameter_list or "DMX" in parameter):

                if (value.find('E') != -1):
                    value = value.replace('E','D')

                if (hasattr(inobj, parameter + 'err')):
                    value_value = getattr(inobj, parameter)
                    error_value = getattr(inobj, parameter + 'err')
                    value = str(value_value)
                    error = ''

                    if (error_value == 0.):
                        
                        error = '0.0'

                    elif (parameter in spin_freq_derivatives):
                        error = str(error_value)

                    else:
                        order = self.order_of_magnitude(getattr(inobj, parameter + 'err'))
                        error = ("{0:." + str(order + 3) + "f}").format(getattr(inobj, parameter + 'err'))

                        if (parameter != 'RAJ' and parameter != 'DECJ'):
                            value = ("{0:." + str(order + 3) + "f}").format(getattr(inobj, parameter))


                    flag = getattr(inobj,parameter+"flag")

                    if (error.find('E') != -1):
                        error = error.replace('E','D')

                    outfile.write("{0:15}    {1:20}  {2:10}    {3:20}\n".format(parameter,value,flag,str(error)))

                elif (hasattr(inobj,parameter+"flag")):
                    flag = getattr(inobj,parameter+"flag")
                    outfile.write("{0:15}    {1:20}  {2:20}\n".format(parameter,value,flag))

                else:
                    outfile.write("{0:15}    {1:20}\n".format(parameter,value))

            else:
                sys.exit("Error: unknown paramater {0}".format(parameter))

        outfile.close()

    def order_of_magnitude(self, error):
        """
        Determines the order of magnitude of the uncertainty.
        """

        log_error = np.log10(error)
        return int(np.fabs(np.round(log_error)))
