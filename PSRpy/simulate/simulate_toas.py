#! /usr/bin/python

from subprocess import Popen, PIPE
from PSRpy.tempo import read_resid2
import numpy as np
import sys

def write_TOAs_to_file(
    toas, 
    toa_uncertainties, 
    frequency_channels,
    n_epochs, 
    n_channels_per_epoch,
    observatory_code = "@",
    output_file="simulated.tim"
):
    """
    Writes simulated TOA data to an ASCII file, assuming Parkes TOA format. 
    """

    fout = open(output_file, "w")
    fout.write("MODE 1\n\n")

    for ii in range(n_epochs):
        for jj in range(n_channels_per_epoch):
            line = " {0:24s} {1:6.4f} {2:20.13f} {3:7.2f} {4:7.2f} {5:>8s}\n".format(
                "fake_data.fits",
                frequency_channels[jj],
                toas[ii],
                0.,
                toa_uncertainties[ii],
                observatory_code
            )

            fout.write(line)

    fout.close()
    
    return 0

def simulate_TOAs(
    parfile,
    bandwidth = 400.,
    central_frequency = 600.,
    epoch_start = 58800., 
    epoch_finish = 58900.,
    jitter_epoch = 5,
    mask_fraction_frequency = 0.1,
    mean_toa_uncertainty = 10.,
    n_epochs = 2,
    n_channels_per_epoch = 1024,
    n_pulses_per_epoch = 1,
    observatory_code = "@",
    output_file = "simulated_toas.tim",
    rms_residual=5.,
    time_range = 365.25,
    use_tempo=True,
    use_tempo2=False
):
    """
    Uses an input parameter file to generate TOAs given a variety of configurable inputs.

    Parameters
    ----------

    parfile : str
        Name of parfile in TEMPO/TEMPO2 format.

    bandwidth : float
        Bandwidth of desired receiver.

    central_frequency : float 
        Central frequency of desired receiver.

    epoch_start : float
        Starting MJD for simulation.

    epoch_finish : float
        Ending MJD for simulation.

    jitter_epoch : int
        The maximum number of days to randomly shift simulated epochs; if non-zero, a 
        random amount of days are added to each simulated in such a way that the quantity
        actual_epoch = original_epoch + numpy.random.uniform(-jitter_epoch, jitter_epoch).

    mask_fraction_frequency : float
        Fraction of channels to randomly zap (i.e., mimic RFI removal)

    mean_toa_uncertainty : float
        Mean value of TOA uncertainty, in microseconds.

    n_epochs : int
        Number of observing epochs to evaluate over the specific time range.

    n_channels_per_epoch : int
        Number of frequency channels across the desired band.

    n_pulses_per_epoch : int
        Number of pulses for which to evaluate TOAs for a given epoch; for each pusle,

    use_tempo : bool 
        If True, use TEMPO for evaluating arrival times.

    use_tempo2 : bool 
        If True, use TEMPO2 for evaluating arrival times.

    Returns
    -------
    """

    n_toas_total = n_epochs * n_pulses_per_epoch * n_channels_per_epoch
    print("Simulating a total of {0} TOAs...".format(n_toas_total))
    print("... number of epochs: {0}".format(n_epochs))
    print("... number of channels per epoch: {0}".format(n_channels_per_epoch))
    print("... number of pulses per epoch: {0}".format(n_channels_per_epoch))

    # first, simulate rough timestamps based on configuration parameters.
    pulse_mjds = np.linspace(epoch_start, epoch_finish, num=n_epochs)
    pulse_mjds += np.random.uniform(-jitter_epoch, jitter_epoch, n_epochs)
    toa_uncertainties = np.fabs(np.random.normal(0., 1., n_epochs)) * rms_residual + mean_toa_uncertainty

    # next, generate the array of frequency channels based on configuration parameters.
    frequency_lower = central_frequency - bandwidth / 2 * (1 - 1 / n_channels_per_epoch)
    frequency_upper = central_frequency + bandwidth / 2 * (1 - 1 / n_channels_per_epoch)
    frequency_channels = np.linspace(frequency_lower, frequency_upper, n_channels_per_epoch)

    # write original, pre-correction TOAs to a file.
    d1 = write_TOAs_to_file(pulse_mjds, toa_uncertainties, frequency_channels, n_epochs, 
             n_channels_per_epoch, observatory_code=observatory_code, 
             output_file="simulated_toas_orig.tim")

    # now, run tempo on these data.
    cmd = ['tempo', '-f', parfile, "simulated_toas_orig.tim"]
    cmd_call = Popen(cmd, stdout=PIPE)
    output, error = cmd_call.communicate()

    for kk in range(3):
    
        # load in output data from initial run.
        toa_data, _ = read_resid2("resid2.tmp")
        corrections = toa_data["residuals"] / 86400.
        #uncertainties = toa_data["toa_uncertainties"]

        # now use the post-fit residuals as corrections, and write a new .tim file.
        pulse_mjds -= corrections
        d1 = write_TOAs_to_file(pulse_mjds, toa_uncertainties, frequency_channels, n_epochs,
                 n_channels_per_epoch, observatory_code=observatory_code, 
                 output_file="simulated_toas_corrected.tim")

        # now, run tempo on these data.
        cmd = ['tempo', '-f', parfile, "simulated_toas_corrected.tim"]
        cmd_call = Popen(cmd, stdout=PIPE)
        output, error = cmd_call.communicate()

    # now, add white noise to corrected data and write to final file.
    pulse_mjds += np.random.normal(0., 1., n_toas_total) * rms_residual * 1e-6 / 86400.

    d1 = write_TOAs_to_file(pulse_mjds, toa_uncertainties, frequency_channels, n_epochs, 
             n_channels_per_epoch, observatory_code=observatory_code, 
             output_file=output_file)


    # clean up.
    cmd = ['rm', 'simulated_toas_orig.tim', 'simulated_toas_corrected.tim']
    cmd_call = Popen(cmd, stdout=PIPE)
    output, error = cmd_call.communicate()

    return 0
