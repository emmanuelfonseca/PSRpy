#! /usr/bin/python -B

import numpy as np
import sys

__all__ = ["Residuals"]

class Residuals:
    """
    A class object that stores and manipulates tempo output. Currently, this module reads in 
    the output of extract.pl (in $TEMPO/util/extract) that contains post-fit residuals and 
    related date in ASCII format. Options are available for plotting data.
    """

    def __init__(self, cfile='c.tmp'):
        """
        Read in tempo output, store as attributes. 
        """

        setattr(self, 'inputfile', cfile)
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

    def average(self, timeavg=1):
        """
        Average residuals collected during an observing session together. Default is to average all 
        TOAs residuals and their uncertainties collected over a 24-hour period.
        """

        return 0


    def plot(self, x='mjd', y='res', reserr=True, info=False, grid=False, 
        resHist=False, savefig=False, figfilename='fig', figfiletype='png', bins=50, 
        fontsize=15, alpha=1, yscalefac=1.5, years=False, plotbothres=True, useclassic=True, ylim=[]):
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
            'year': 'Year',
            'res': r'Post-fit Residual ($\mu$s)',
            'res_P': r'Post-fit Residual (fractional spin-period)',
            'orb_phase': 'Orbital Phase',
            'uncertainty': "Post-fit Residual Uncertainty ($\mu$s)"
        }

        # if desired, read in info.tmp file.
        info_flags = []
        if (info):
            try:
                info_flags = np.loadtxt(info, dtype=str)
            except:
                sys.exit("Cannot read info-flag file.")

        # if desired, use classic look for matplotlib.
        if useclassic:
            plt.style.use('classic')

        fig, ax = plt.subplots()

        # generate the desired plot.
        if (resHist):
            ax.hist(getattr(self, 'res'), bins, alpha=alpha)
        else:

            x_data = getattr(self, x)
            y_data = getattr(self, y)
            # if desired x-axis is time in years, convert.
            if (x == 'mjd' and years):
                x = 'year'
                x_data = (x_data - 53005.) / 365.25 + 2004.
            if (y == 'res' and reserr):
                yerr_data = self.uncertainty
                if (info):
                    for label in np.unique(info_flags):
                        x_data_int = x_data[(np.where(info_flags == label))[0]]
                        y_data_int = y_data[(np.where(info_flags == label))[0]]
                        yerr_data_int = yerr_data[(np.where(info_flags == label))[0]]
                        ax.errorbar(x_data_int, y_data_int, yerr=yerr_data_int, fmt='+')
                        
                else:
                    ax.errorbar(x_data, y_data, yerr=yerr_data, fmt='+')
            else:
                plt.plot(x_data, y_data, 'b+')

        ax.set_xlabel(axlabel[x], fontproperties=font, fontsize=fontsize)
        ax.set_ylabel(axlabel[y], fontproperties=font, fontsize=fontsize)

        # add grids, if desired.
        if (grid): 
            ax.grid()

        # now, set y-axis limits. 
        ax.set_ylim(np.min(y_data) * yscalefac, np.max(y_data) * yscalefac)

        if (len(ylim) == 2):
            ax.ylim(ylim)

        if plotbothres:
            ax2 = ax.twinx()
            ax2.set_ylabel(axlabel['res_P'], fontproperties=font, fontsize=fontsize)
            ax2.set_ylim(np.min(self.res_P) * yscalefac, np.max(self.res_P) * yscalefac)
            plt.tight_layout()

        # save figure in png format, if desired.
        if savefig:
            plt.savefig(figfilename + '.' + figfiletype, fmt=figfiletype)

        plt.show()

    def stats(self, verbose=False):
        """
        Compute statistics for timing residuals.
        """

        weight = 1 / (self.uncertainty)**2 / np.sum(1 / (self.uncertainty)**2)

        setattr(self, 'Ntoa', len(self.res))
        setattr(self, 'rms', np.sqrt(np.sum((self.res)**2) / self.Ntoa))
        setattr(self, 'rms_weighted', np.sqrt(np.sum(weight * (self.res)**2)))
        setattr(self, 'median_uncertainty', np.median(self.uncertainty))

        if (verbose):
            print "Residual statistics for {0}".format(self.inputfile)
            print "    * number of residuals: {0}".format(self.Ntoa)
            print "    * RMS residual: {0:.3f} microsec".format(self.rms)
            print "    * ... weighted: {0:.3f} microsec".format(self.rms_weighted)
