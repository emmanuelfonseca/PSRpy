#! /usr/bin/python

import matplotlib.pyplot as plt
import numpy.fft as nf
import numpy as np
import sys

def tay92_equation_A7(tau, amps, ampd, angs, angd, k):
    """
    A function that returns the value of equation A7 in Taylor (1992).
    """

    return np.sum(k * amps * ampd * np.sin(-angs + angd + k * (tau * 2 * np.pi)))

def tay92_equation_A9(tau, amps, ampd, angs, angd, k):
    """
    A function that returns the scale factor between the template and data profiles.
    See Equation 9 of Taylor (1992).
    """

    return np.sum(amps * ampd * np.cos(-angs + angd + k * (tau * 2 * np.pi))) / np.sum(amps**2)

def fftfit(data_prof, tmp_prof, tolerance=1e-12):
    """
    A function that uses the Taylor (1992) method to determine the shift between a template 
    profile with a noisey profile.    
    """

    tmp_fft = nf.rfft(tmp_prof)
    amp_tmp = np.absolute(tmp_fft)
    phs_tmp = np.angle(tmp_fft)
    k_tmp = np.linspace(0, len(amp_tmp)-1, num=len(amp_tmp))

    dat_fft = nf.rfft(data_prof)
    amp_dat = np.absolute(dat_fft)
    phs_dat = np.angle(dat_fft)
    k_dat = np.linspace(0, len(amp_dat)-1, num=len(amp_dat))

    # compute Equation A7 in Taylor (1992) for different values of index k.

    phs_array = np.linspace(0, 1, num=len(k_dat))
    A7vals = np.zeros(len(phs_array))
    count = 0

    for phs in phs_array:
        A7vals[count] = tay92_equation_A7(phs, amp_tmp[1:], amp_dat[1:], phs_tmp[1:], phs_dat[1:], k_tmp[1:])
        count += 1

    # implement Brent's method for finding the zero to Equation A7.

    c, f_best = 0, 0
    a = (phs_array[np.where(A7vals == max(A7vals))])[0]
    b = (phs_array[np.where(A7vals == min(A7vals))])[0]

    while (np.fabs(b - a) > tolerance):
        c = (a + b) / 2
        f_best = tay92_equation_A7(c, amp_tmp[1:], amp_dat[1:], phs_tmp[1:], phs_dat[1:], k_tmp[1:])

        if (f_best < 0.):
            b = c
        else:
            a = c

    best_shift = c
    best_scale = tay92_equation_A9(best_shift, amp_tmp[1:], amp_dat[1:], phs_tmp[1:], phs_dat[1:], k_tmp[1:])

    # compute relative offset of profile baselines.
    best_offset = (amp_dat[0] - best_scale * amp_tmp[0]) / len(data_prof)

    return best_offset, best_shift, best_scale

def fftshift(prof, tau=0., scale=1.0, baseline=0.):
    """
    A function that shifts an input profile by an amount tau in the Fourier domain.
    """

    ftp = nf.rfft(prof)
    ampf = np.absolute(ftp)
    angf = np.angle(ftp)
    kf   = np.linspace(0., len(ftp)-1., num=len(ftp))

    return np.fft.irfft(ampf * np.exp (1j * (angf + kf * (tau * 2 * np.pi))) * scale) + baseline

def diffprof(data_prof, tmp_prof):
    """
    A function that computes the difference between two profiles; the difference is taken 
    in the Fourier domain.
    """

    baseline, shift, fac = fftfit(data_prof, tmp_prof)
    tmp_prof_shifted = fftshift(tmp_prof, tau=-shift, scale=fac, baseline=baseline)
    
    return data_prof - tmp_prof_shifted
