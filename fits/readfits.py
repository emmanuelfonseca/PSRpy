#! /usr/bin/python

from PSRpy.dm import dm_time_delay
from astropy.io import fits
import PSRpy.fft as ft
import matplotlib.pyplot as plt
import numpy as np
import sys

class ReadFits():
    """
    Defines a class that stores data and header info from an input PSRFITS file.
    """

    def __init__(self, infile):
        """
        Reads in FITS file using the astropy.io.fits module.

        Input:
            - FITS file that uses PSRFITS format.
        """
        hdulist = fits.open(infile)

        # store info from main header of fits file.
        date, starttime = (hdulist[0].header['DATE-OBS']).split('T')
        setattr(self, 'inputfile', infile)
        setattr(self, 'epoch', date)
        setattr(self, 'starttime', starttime)
        setattr(self, 'projectID', hdulist[0].header['PROJID'])
        setattr(self, 'telescope', hdulist[0].header['TELESCOP'])
        setattr(self, 'observer', hdulist[0].header['OBSERVER'])
        setattr(self, 'source', hdulist[0].header['SRC_NAME'])
        setattr(self, 'receiver', hdulist[0].header['FRONTEND'])
        setattr(self, 'backend', hdulist[0].header['BACKEND'])
        setattr(self, 'mode', hdulist[0].header['OBS_MODE'])

        # store key numbers from 'history' binary table extension.
        setattr(self, 't_bin', (hdulist[1].data['TBIN'])[0])
        dedisp = (hdulist[1].data['DEDISP'])[0]
        if dedisp:
            setattr(self, 'dedisp', True)
        else:
            setattr(self, 'dedisp', False)

        # store key header info, pulsar data from 'subint data' table extension.
        data  = hdulist[4].data['DATA']
        setattr(self, 'data', data)
        setattr(self, 'channel_freqs', np.array((hdulist[4].data['DAT_FREQ'])[0, :]))
        setattr(self, 'n_ints', len(data[:, 0, 0 ,0]))
        setattr(self, 'n_bins', hdulist[4].header['NBIN'])
        setattr(self, 'n_chan', hdulist[4].header['NCHAN'])
        setattr(self, 'n_pol', hdulist[4].header['NPOL'])
        setattr(self, 'dm', hdulist[4].header['DM'])

    def info(self):
        """
        Prints basic info of input file.
        """

        print "Basic header info for file '{0}':".format(self.inputfile)
        print "    * Date of observation: {0} (started at {1})".format(self.epoch, self.starttime)
        print "    * Project ID and telescope: {0} ({1})".format(self.projectID, self.telescope)
        print "    * Observer: {0}".format(self.observer)
        print "    * Source: {0}".format(self.source)
        print "    * Receiver: {0}".format(self.receiver)
        print "    * Backend: {0}".format(self.backend)
        print "    * Mode: {0}".format(self.mode)

    def heatmap_phase_frequency(self, pol=0, reference_freq=430., ignore_chans=[], ignore_subints=[]):
        """
        Computes full-sum heat map in frequency and orbital phase.

        Inputs:
            - pol            : polarization proflie (pol=0 -> total intensity)
            - reference_freq : reference frequency for dispersion removal (MHz)
            - ignore_chans   : channels to ignore (python list of integers)
            - ignore_subints : subints to ignore (python list of integers)

        Outputs:
            - freq_good           : array of frequencies that excluded the ignored ones.
            - phase_frequency_map : NumPy array of signal, with shape (n_freq_good x n_bins)
        """

        period_topo = self.t_bin * self.n_bins
        phase_frequency_map = np.zeros((self.n_chan - len(ignore_chans), self.n_bins))
        freq_good = np.zeros(self.n_chan - len(ignore_chans))
        count = 0

        for ii in range(self.n_chan):
            curr_prof = np.zeros(self.n_bins)
            if (ii not in ignore_chans):
                for jj in range(self.n_ints):
                    curr_prof += self.data[jj, pol, ii, :] / self.n_ints
                if (not self.dedisp):
                    shift_dm = dm_time_delay(self.dm, self.channel_freqs[ii], reference_freq) / period_topo
                    phase_frequency_map[count, :] = ft.fftshift(curr_prof, tau=shift_dm)
                else:
                    phase_frequency_map[count, :] = curr_prof
                freq_good[count] = self.channel_freqs[ii]
                count += 1

        return freq_good, phase_frequency_map

    def heatmap_phase_time(self, pol=0, reference_freq=430., ignore_chans=[], ignore_subints=[]):
        """
        Computes full-sum heat map in time and pulse phase.

        Inputs:
            - pol            : polarization proflie (pol=0 -> total intensity)
            - reference_freq : reference frequency for dispersion removal (MHz)
            - ignore_chans   : channels to ignore (python list of integers)
            - ignore_subints : subints to ignore (python list of integers)

        Outputs:
            - phase_time_map : NumPy array of signal, with shape (n_ints x n_bins)
        """

        phase_time_map = np.zeros((n_ints, n_bins))
        period_topo = self.t_bin * self.n_bins

        for kk in range(self.n_ints):
            n_chan_good = 0
            for ll in range(self.n_chan):
                if (ll not in ignore_chans):
                    # if not done so already, de-disperse signal.
                    if (not self.dedisp):
                        shift_dm = dm_time_delay(self.dm, freq[ll], reference_freq) / period_topo
                        phase_time_map[kk, :] += ft.fftshift(self.data[kk, 0, ll, :], tau=shift_dm) 
                    else:
                        phase_time_map[kk, :] += self.data[kk, 0, ll, :]
                    n_chan_good += 1
            phase_time_map[kk, :] /= n_chan_good

        return phase_time_map

    def full_summed_profile(self, pol=0, reference_freq=430., ignore_chans=[], ignore_subints=[]):
        """
        Computes full-sum heat map in time and pulse phase.

        Inputs:
            - pol            : polarization proflie (pol=0 -> total intensity)
            - reference_freq : reference frequency for dispersion removal (MHz)
            - ignore_chans   : channels to ignore (python list of integers)
            - ignore_subints : subints to ignore (python list of integers)

        Outputs:
            - phase_time_map : NumPy array of signal, with shape (n_ints x n_bins)
        """

        full_sum_profile = np.zeros(self.n_bins)
        freqs, heat_map = self.heatmap_phase_frequency(pol=pol, reference_freq=reference_freq, 
                                                       ignore_chans=ignore_chans, ignore_subints=ignore_subints)

        for kk in range(len(freqs)):
            full_sum_profile += heat_map[kk, :]

        return full_sum_profile
