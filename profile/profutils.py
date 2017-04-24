#! /usr/bin/python

import numpy as np
#import pywt

#def smooth_wavelet1D(indat, wave='sym8', nlevel=5, ncycle=10, threshtype='hard'):
#    """
#    Compute the wavelet-denoised version of an input profile/waveform. 
#    
#    Required arguments: 
#        - 'indat' = input data in 1-dimensional array. 
#
#    Default arguments: 
#        - 'wave' = mother wavelet (default 'sym8') 
#        - 'nlevel' = number of decomposition levels (default '5') 
#        - 'ncycle' = number of circulant averages to compute (default '10') 
#        - 'threshtype' = type of wavelet thresholding (default 'hard')
#    """

#    nbins = len(indat)
#    data  = np.zeros(nbins)

#    # carry out translation-invariant wavelet denoising.
#    for j in range(ncycle):
#        m = j - ncycle/2 - 1
#        coeffs = pywt.wavedec(np.roll(indat,m),wave,level=nlevel)
#        # get threshold value.
#        lopt = np.median(np.fabs(coeffs[1]))/0.6745*np.sqrt(2.*np.log(nbins))
#        # now do wavelet thresholding.
#        for k in range(1,nlevel+1):
#        # hard threshold.
#            if (threshtype == 'hard'):
#                (coeffs[k])[np.where(np.fabs(coeffs[k]) < lopt)] = 0.
#        # soft threshold. 
#        else:
#            (coeffs[k])[np.where(np.fabs(coeffs[k]) < lopt)] = 0.
#            (coeffs[k])[np.where(coeffs[k] > lopt)] = (coeffs[k])[np.where(coeffs[k] > lopt)]+lopt
#            (coeffs[k])[np.where(coeffs[k] < -lopt)] = (coeffs[k])[np.where(coeffs[k] < -lopt)]-lopt
#        # reconstruct data.
#        data = data + np.roll(pywt.waverec(coeffs,wave),-m)   
  
#    return data/np.float(ncycle)

#def smooth_wavelet2D(indat, wave='sym8', nlevel=5, ncycle=10, threshtype='hard'):
#    """
#    Compute the wavelet-denoised version of a set of profiles/waveforms.
#  
#    Required argument:
#        - 'indat' = (nbins x nchans) data array.
#  
#    Default arguments: 
#        - 'wave' = mother wavelet (default 'sym8') 
#        - 'nlevel' = number of decomposition levels (default '5') 
#        - 'ncycle' = number of circulant averages to compute (default '10') 
#        - 'threshtype' = type of wavelet thresholding (default 'hard')
#    """

#    nbins = len(indat[:,0])
#    nchan = len(indat[0,:])

#    outdat = np.zeros((nbins,nchan))
  
#    # smooth each channel.
#    for i in range(nchan):
#        prof = indat[:,i]
#        data = 0.
#        # carry out translation-invariant wavelet denoising.
#        for j in range(ncycle):
#            m = j - ncycle/2 - 1
#            coeffs = pywt.wavedec(np.roll(prof,m),wave,level=nlevel)
#            # get threshold value.
#            lopt = np.median(np.fabs(coeffs[1]))/0.6745*np.sqrt(2.*np.log(nbins))
#            # now do wavelet thresholding.
#            for k in range(1,nlevel+1):
#                # hard threshold.
#                if (threshtype == 'hard'):
#                    (coeffs[k])[np.where(np.fabs(coeffs[k]) < lopt)] = 0.
#                # or soft threshold.
#                else:
#                    (coeffs[k])[np.where(np.fabs(coeffs[k]) < lopt)] = 0.
#                    (coeffs[k])[np.where(coeffs[k] > lopt)] = (coeffs[k])[np.where(coeffs[k] > lopt)]+lopt
#                    (coeffs[k])[np.where(coeffs[k] < -lopt)] = (coeffs[k])[np.where(coeffs[k] < -lopt)]-lopt
#            # reconstruct data.
#            data = data + np.roll(pywt.waverec(coeffs,wave),-m)
#        # save averaged profile.
#        outdat[:,i] = data/np.float(ncycle)
#
#    # return smoothed portrait.
#    return outdat

def smooth_pca(indat, ncomp=5):
    """
    Compute the PCA-denoised version of a set of profiles/waveforms.
 
    Required argument:
        - indat = (nbins x nchan) data array.

    Default argument:
        - ncomp = number of principal components to use (default '5')
    """
  
    nmes = len(indat[0,:])
    ndim = len(indat[:,0])

    print 'PCA on data with {0} dimensions, {1} measurements:'.format(ndim,nmes)

    # subtract mean from each set of measurements.
    dmean  = np.mean(indat,axis=1)
    matrix = (indat.T-dmean).T

    # compute covariance matrix.
    cov = np.cov(indat)

    # compute eigenvalues/vectors of cov, and order them.
    eigval, eigvec = np.linalg.eigh(cov)
    ind = (np.argsort(eigval))[::-1]
    eigval, eigvec = eigval[ind], eigvec[:,ind]
 
    # transform to new, smoothed data using desired number of components.
    eigvec = eigvec[:,range(ncomp)] 
    finaldata = (np.dot(eigvec,np.dot(eigvec.T,matrix)).T+dmean).T
  
    return finaldata
