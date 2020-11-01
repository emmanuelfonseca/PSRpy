from ..const import k_DM
import astropy.units as u
import numpy as np

def dm_delay(dm, frequency1, frequency2=np.infty):
    """
    Computes the time delay between frequency channels due to electromagnetic 
    dispersion, given a diserpsion measure.

    Parameters
    ----------
    dm : float
        dispersion measure, in units of pc cm^{-3}

    frequency1 : float 
        central frequency of a channel for which to evaluate delay, in units of MHz

    frequency2 : float, optional
        central freqneyc of a different channel, in units of MHz (default: infinity)

    Returns
    -------
    delay : 
    """

    dm_in = dm * u.pc / u.cm / u.cm / u.cm    
    freq1_in = frequency1 * u.MHz
    freq2_in = frequency2 * u.MHz
    delay = k_DM * dm_in * (1 / freq1_in**2 - 1 / freq2_in**2)

    return delay

def dm_smearing_delay(dm, frequency, frequency_width):
    """
    Computes the time delay associated with smearing of a dispersed pulse resolved into 
    finite frequency bins. 

    Parameters
    ----------

    """

    dm_in = dm * u.pc / u.cm / u.cm / u.cm    
    freq_in = frequency * u.MHz
    freq_width_in = frequency_width * u.MHz
    delay = 2 * k_DM * dm_in * freq_width_in / freq_in**3

    return delay
