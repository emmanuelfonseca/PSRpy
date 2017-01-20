import numpy as np

def dm_time_delay(dm, freq1, freq2):
    """
    Computes the time delay between frequency channels due to electromagnetic 
    dispersion, given a diserpsion measure.

    Inputs:
        - dm    : dispersion measure (cm^{-3} pc)
        - freq1 : central frequency of a channel (MHz)
        - freq2 : central freqneyc of a different channel (MHz)
    """

    constant = 4.148808e3
    return dm * constant * (1 / freq1**2 - 1 / freq2**2)
