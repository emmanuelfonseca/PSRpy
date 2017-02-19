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
        date, starttime = (hdulist['PRIMARY'].header['DATE-OBS']).split('T')
        setattr(self, 'inputfile', infile)
        setattr(self, 'epoch', date)
        setattr(self, 'starttime', starttime)
        setattr(self, 'projectID', hdulist['PRIMARY'].header['PROJID'])
        setattr(self, 'telescope', hdulist['PRIMARY'].header['TELESCOP'])
        setattr(self, 'observer', hdulist['PRIMARY'].header['OBSERVER'])
        setattr(self, 'source', hdulist['PRIMARY'].header['SRC_NAME'])
        setattr(self, 'receiver', hdulist['PRIMARY'].header['FRONTEND'])
        setattr(self, 'backend', hdulist['PRIMARY'].header['BACKEND'])
        setattr(self, 'mode', hdulist['PRIMARY'].header['OBS_MODE'])

        # store key numbers from 'history' binary table extension.
        setattr(self, 't_bin', (hdulist['HISTORY'].data['TBIN'])[0])
        setattr(self, 'dedisp', hdulist['HISTORY'].data['DEDISP'])

        # store key header info, pulsar data from 'subint data' table extension.
        data  = hdulist['SUBINT'].data['DATA']
        scale = hdulist['SUBINT'].data['DAT_SCL']
        offset = hdulist['SUBINT'].data['DAT_OFFS']
        setattr(self, 'scale', np.array(scale, dtype=np.float32))
        setattr(self, 'offset', np.array(offset, dtype=np.float32))
        setattr(self, 'channel_freqs', np.array((hdulist['SUBINT'].data['DAT_FREQ'])[0, :]))
        setattr(self, 'n_ints', len(data[:, 0, 0 ,0]))
        setattr(self, 'n_bins', hdulist['SUBINT'].header['NBIN'])
        setattr(self, 'n_bits', hdulist['SUBINT'].header['NBITS'])
        setattr(self, 'n_chan', hdulist['SUBINT'].header['NCHAN'])
        setattr(self, 'n_pol', hdulist['SUBINT'].header['NPOL'])
        setattr(self, 'signint', hdulist['SUBINT'].header['SIGNINT'])
        setattr(self, 'dm', hdulist['SUBINT'].header['DM'])

        setattr(self, 'data', np.array(data, dtype=np.float))

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

    def dedisperse(self, reference_freq=430.):
        """
        De-disperses fold-mode data in given fits file.
        """

        period_topo = self.t_bin * self.n_bins
        shift_dm = 0.

        # loop over polarization, channel and subint.

        for pol in range(self.n_pol):
            freq_idx = 0
            for freq in self.channel_freqs:
                # compute DM time delay.
                shift_dm = dm_time_delay(self.dm, freq, reference_freq) / period_topo
                for subint in range(self.n_ints):
                    # dedisperse by shifting to common phase.
                    self.data[subint, pol, freq_idx, :] = ft.fftshift(self.data[subint, pol, freq_idx, :], tau=shift_dm)
                freq_idx += 1

    def downsample(self, freqs, chan_data, new_bins, pol=0):
        """
        Downsamples each sub-integration to a new number of bins.

        TODO: re-write this so that it directly edits the fits object.
        """

        factor = self.n_bins / new_bins

        if (float(factor).is_integer()):
            chan_data_downsampled = np.zeros((len(freqs), new_bins))      

            for freq in range(len(freqs)):
                cbin = 0
                for pbin in range(new_bins):
                    chan_data_downsampled[freq, pbin] = np.mean(chan_data[freq, cbin:(pbin+1)*factor])
                    cbin += factor

            return chan_data_downsampled

        else:
        
            print "WARNING: requested bin number not even number! Proceeding with original data..."
            return chan_data

    def remove_baseline(self, phase_range=[0.7, 0.8]):
        """
        Removes non-zero baseline, such that off-pulse RMS intensity is zero.
        """
        
        pulse_phase = np.linspace(0, self.n_bins-1, num=self.n_bins) / np.float(self.n_bins)
        bl_idx = np.where(np.logical_and(pulse_phase >= phase_range[0], pulse_phase <= phase_range[1]))[0]

        for pol in range(self.n_pol):
            for freq in range(self.n_chan):
                for subint in range(self.n_ints):
                    self.data[subint, pol, freq, :] -= np.mean(self.data[subint, pol, freq, bl_idx])

    def rescale(self):
        """
        Re-scales data given the scale and offset parameters.
        """
        
        for pol in range(self.n_pol):
            for freq in range(self.n_chan):
                for subint in range(self.n_ints):
                    self.data[subint, pol, freq, :] *= self.scale[subint, pol * self.n_chan + freq]
                    self.data[subint, pol, freq, :] += self.offset[subint, pol * self.n_chan + freq]

    def shift_phase(self, shift):
        """
        Shifts all profiles in pulse phase.
        """

        for pol in range(self.n_pol):
            for freq in range(self.n_chan):
                for subint in range(self.n_ints):
                    self.data[subint, pol, freq, :] = ft.fftshift(self.data[subint, pol, freq, :], tau=shift)

    def heatmap_phase_frequency(self, pol=0, reference_freq=430., shift_by=0.2, ignore_chans=[], 
                                ignore_subints=[], dedisp=False, rescale=False, rm_baseline=False, shift=False):
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

        phase_frequency_map = np.zeros((self.n_chan - len(ignore_chans), self.n_bins))
        freq_good = np.zeros(self.n_chan - len(ignore_chans))
        count = 0

        if dedisp:
            self.dedisperse(reference_freq=reference_freq)

        if shift:
            self.shift_phase(shift_by)

        if rescale:
            self.rescale()

        if rm_baseline:
            self.remove_baseline()

        for ii in range(self.n_chan):
            curr_prof = np.zeros(self.n_bins)
            if (ii not in ignore_chans):
                for jj in range(self.n_ints):
                    curr_prof = np.array(self.data[jj, pol, ii, :]) / self.n_ints
                    phase_frequency_map[count, :] += curr_prof
                freq_good[count] = self.channel_freqs[ii]
                count += 1

        return freq_good, phase_frequency_map

    def heatmap_phase_time(self, pol=0, reference_freq=430., ignore_chans=[], ignore_subints=[], 
                           dedisp=False, rescale=False, rm_baseline=False, shift=False):
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

        phase_time_map = np.zeros((self.n_ints, self.n_bins))
        period_topo = self.t_bin * self.n_bins

        if dedisp:
            self.dedisperse(reference_freq=reference_freq)

        if rescale:
            self.rescale()

        if rm_baseline:
            self.remove_baseline()

        for kk in range(self.n_ints):
            n_chan_good = 0
            for ll in range(self.n_chan):
                if (ll not in ignore_chans):
                    phase_time_map[kk, :] += self.data[kk, 0, ll, :]
                    n_chan_good += 1
            phase_time_map[kk, :] /= n_chan_good

        return phase_time_map

    def full_summed_profiles(self, reference_freq=430., shift_by=0.2, ignore_chans=[], ignore_subints=[], 
                             dedisp=False, rm_baseline=False, rescale=False, shift=False):
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

        full_sum_profiles = np.zeros((self.n_bins, self.n_pol))

        if dedisp:
            self.dedisperse()

        if shift:
            self.shift_phase(shift_by)

        if rescale:
            self.rescale()

        if rm_baseline:
            self.remove_baseline()

        for pol in range(self.n_pol):
            for freq in range(self.n_chan):
                for subint in range(self.n_ints):
                    full_sum_profiles[:, pol] += self.data[subint, pol, freq, :] 


        return full_sum_profiles / np.float(self.n_chan * self.n_ints)
