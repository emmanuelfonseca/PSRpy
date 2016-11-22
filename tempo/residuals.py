#! /usr/bin/python -B

import numpy as np
import sys

class Residuals():
    """
    A class object that stores and manipulates tempo output. Currently, this module reads in 
    the output of extract.pl (in $TEMPO/util/extract) that contains post-fit residuals and 
    related date in ASCII format. Options are available for plotting data.
    """

    def __init__(self, cfile='c.tmp'):
        """
        Read in tempo output, store as attributes.
        """

        data = np.loadtxt(cfile)
        setattr(self, 'daynum', data[:, 0])
        setattr(self, 'res_P', data[:, 1])
        setattr(self, 'res', data[:, 2])
        setattr(self, 'uncertainty', data[:, 3])
        setattr(self, 'res_prefit', data[:, 4])
        setattr(self, 'orb_phase', data[:, 5])
        setattr(self, 'mjd', data[:, 6])
        setattr(self, 'ptype', data[:, 7])
        setattr(self, 'serial', data[:, 8])
        setattr(self, 'daynum_year', data[:, 9])
        setattr(self, 'freq', data[:, 10])
        setattr(self, 'timeofday', data[:, 11])
        setattr(self, 'serial2', data[:, 12])
        setattr(self, 'mjd2', data[:, 13])
        setattr(self, 'weight', data[:, 14])

    def plot(self, x='mjd', y='res', reserr=True, info=False, grid=False, 
        resHist=False, savefig=False, figfilename='fig', figfiletype='png', bins=50, 
        fontsize=15, alpha=1):
        """
        Plot data of choice. 
        """

        from matplotlib.font_manager import FontProperties
        import matplotlib.pyplot as plt
        font = FontProperties()
        font.set_name('serif')

        # define axis-label dictionary.
        axlabel = {
            'mjd': 'MJD',
            'res': r'Post-fit Residual ($\mu$s)',
            'orb_phase': 'Orbital Phase'
        }

        # if desired, read in info.tmp file.
        info_flags = []
        if (info):
            try:
                info_flags = np.loadtxt(info, dtype=str)
            except:
                sys.exit("Cannot read info-flag file.")

        # generate the desired plot.
        if (resHist):
            plt.hist(getattr(self, 'res'), bins, alpha=alpha)
        else:
            x_data = getattr(self, x)
            y_data = getattr(self, y)
            if (y == 'res' and reserr):
                yerr_data = self.uncertainty
                if (info):
                    for label in np.unique(info_flags):
                        x_data_int = x_data[(np.where(info_flags == label))[0]]
                        y_data_int = y_data[(np.where(info_flags == label))[0]]
                        yerr_data_int = yerr_data[(np.where(info_flags == label))[0]]
                        plt.errorbar(x_data_int, y_data_int, yerr=yerr_data_int, fmt='+')
                else:
                    plt.errorbar(x_data, y_faya, yerr=yerr_data, fmt='+')
            elif(y == 'res'):
                plt.plot(x_data, y_data)
        plt.xlabel(axlabel[x], fontproperties=font, fontsize=fontsize)
        plt.ylabel(axlabel[y], fontproperties=font, fontsize=fontsize)

        # add grids, if desired.
        if (grid): 
            plt.grid()
        
        # save figure in png format, if desired.
        if (savefig):
            plt.savefig(figfilename + '.' + figfiletype, fmt=figfiletype)

        plt.show()
