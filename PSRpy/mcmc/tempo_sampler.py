#! /usr/bin/python

from PSRpy.parfile import Parfile, PrintPar
from astropy.coordinates import Angle
from subprocess import Popen, PIPE
from copy import deepcopy
import astropy.units as u
import numpy as np
import sys
import re

class TempoMCMC():

    def __init__(self, parfile, timfile, marginalize=[]):
        """
        Initialize sampler class.
        """

        self.parfile = parfile
        self.timfile = timfile

        self.par_object = Parfile(self.parfile)

        current_parameters = self.par_object.fit_parameters
        parameters_for_mcmc = []

        for curr_par in  current_parameters:
            if (curr_par not in marginalize):
                parameters_for_mcmc.append(curr_par)

        self.parameters = parameters_for_mcmc
        self.n_params = len(self.parameters)
        self.priors = {}
        self.posteriors = {}
        self.step_sizes = {}
        self.tempomodel = []

    def log_likelihood(self, parfile, tempo=True, gls=False, pulse_numbers=None):
        """
        Computes the likelihood.

        Inputs
        ------

        parfile : str
            Pulsar-timing parameter file, in TEMPO format.

        tempo : bool
            If set to True, use TEMPO to compute chi-squared value; otherwise, use TEMPO2.

        gls : bool 
            If set to True, use generalized least-squares (GLS) fitting in TEMPO.

        pulse_numbers : str
            ASCII file that contains pulse numbers; optional.
        """

        tempo_command = []

        if tempo:

            if gls:

                if (pulse_numbers is not None):
                    tempo_command = ['tempo', '-G', '-ni', pulse_numbers, '-f ', parfile, self.timfile]

                else:
                    tempo_command = ['tempo', '-G', '-f ', parfile, self.timfile]

            else:

                if (pulse_numbers is not None):
                    tempo_command = ['tempo', '-ni', pulse_numbers, '-f ', parfile, self.timfile]

                else:
                    tempo_command = ['tempo', '-f ', parfile, self.timfile]


        temporun = Popen(tempo_command, stdout=PIPE)
        output, error = temporun.communicate()

        m = re.search(r'Chisqr/nfree\W+(\d+\.\d+)\W+(\d+)\W+(\d+\.\d+)', output)
        chisq = np.float(m.group(3)) * np.float(m.group(2))

        return -0.5 * chisq

    def plot(self, n_bins=100):
        """
        Plots MCMC posterior data.
        """

        from matplotlib.font_manager import FontProperties
        import matplotlib.pyplot as plt

        font = FontProperties()
        font.set_name('serif')

        # plot chains versus iterations.

        for num, parameter in zip(range(self.n_params), self.parameters):
            #3plt.subplot(self.n_params, 1, num + 1)
            if (parameter == 'F0' or parameter == 'PB' or parameter == 'T0'):
                plt.plot(self.posteriors[parameter] - np.mean(self.posteriors[parameter]), 'b-')
 
            else:
                plt.plot(self.posteriors[parameter], 'b-')
            plt.savefig('chain.' + parameter + '.png', fmt='png')
            plt.gcf().clear()

        #plt.show()
        ##sys.exit()

        #plt.figure(2)

        #from matplotlib.gridspec import GridSpec
        #gs = GridSpec(self.n_params, self.n_params)

        #count = 0

        #for num1, parameter1 in zip(range(self.n_params), self.parameters):
        #    for num2, parameter2 in zip(range(self.n_params), self.parameters):
        #        if (num1 == num2):
        #            plt.subplot(gs[num1, num2])
        #            plt.hist(self.posteriors[parameter1], bins=n_bins, histtype='step')
        #            plt.xticks([])
        #            plt.yticks([])
        #        elif (num1 > num2):
        #            plt.subplot(gs[num1, num2])
        #            hist2D, xedges, yedges = np.histogram2d(self.posteriors[parameter2], self.posteriors[parameter1], bins=50)
        #            dx = xedges[1] - xedges[0]
        #            dy = yedges[1] - yedges[0]
        #            x = np.linspace(xedges[0] + dx, xedges[len(xedges)-2] + dx, num=len(xedges)-1)
        #            y = np.linspace(yedges[0] + dy, yedges[len(yedges)-2] + dy, num=len(yedges)-1)
        #            plt.pcolormesh(x, y, hist2D, cmap='Blues')
        #        else:
        #            pass
        #        count += 1

        #plt.show()

    def prior(self, prior_type, all_priors=True, parameter_name=None):
        """
        Sets type of prior distribution for one or all parameters.

        Inputs
        ------

        prior_type : str
            Name of prior to use for one or all parameters.

        all_priors : bool
            If true, all priors will be set to 'prior_type'.

        parameter_name : str
            Name of parameter that will have prior with type 'prior_type'.
        """

        if all_priors:
        
            for parameter in self.parameters:
                self.priors[parameter] = prior_type

        else:
        
            if (parameter_name not in self.parameters):
                sys.exit("Error: option parameter_name not in list of allowed parameters!")
            
            self.priors[parameter_name] = prior_type

    def run(self, burn=500, efac=1, n_iterations=1000, thin=0):
        """
        Runs the MCMC sampler.

        Inputs
        ------

        burn : int
            Number of iterations to ignore (i.e. burn). Default is 500.

        n_iterations : int
            Number of iterations for sampler. Default is 1000.

        thin : int
            Thinning factor. Default is 0.

        TODO: re-write this for TEMPO/parfile; implement thinning.
        """

        samples = np.zeros((n_iterations, self.n_params))
        raj_samples = []
        decj_samples = []
        success_count = 0
        self.par_object.write(outfile='mcmc_current.par')

        for iternum in range(n_iterations):

            current_par = Parfile('mcmc_current.par')
            setattr(self, 'par_object', deepcopy(current_par))
            self.step(factor=efac)
            self.par_object.write(outfile='mcmc_proposal.par')
            ll_current = self.log_likelihood('mcmc_current.par')
            ll_proposal = self.log_likelihood('mcmc_proposal.par')
            log_prior = np.log(np.random.uniform(0., 1.))

            if (log_prior < np.min([0., ll_proposal - ll_current])):

                mv_command = Popen(['mv', 'mcmc_proposal.par', 'mcmc_current.par'], stdout=PIPE)
                run_mv, out_mc = mv_command.communicate()
                success_count += 1
                print(success_count, iternum + 1)

            else:
                setattr(self, 'par_object', deepcopy(current_par))

            count = 0
            for parameter in self.parameters:

                if (parameter == 'RAJ' and parameter != 'DECJ'):
                    ra = Angle(getattr(self.par_object, 'RAJ'), unit=u.hour)
                    samples[iternum, count] = ra.deg
                elif (parameter == 'DECJ'):
                    dec = Angle(getattr(self.par_object, 'DECJ'), unit=u.deg)
                    samples[iternum, count] = dec.deg
                else:
                    samples[iternum, count] = getattr(self.par_object, parameter)

                count += 1

        self.accepted = success_count
        
        for num, parameter in zip(range(self.n_params), self.parameters):
                self.posteriors[parameter] = samples[burn:, num]

    def step(self, factor=1):
        """
        Steps a set of parameters, given information on prior.

        Inputs
        ------

        factor : int
            Multiplicative value applied to every uncertainty.

        TODO: adapt this module for Parfile step/write functionality.
        """

        n_dim = len(self.parameters) / 2.4**2

        for parameter in self.parameters:
            value = getattr(self.par_object, parameter)
            err = self.step_sizes[parameter] * factor / n_dim

            if (err != 0.):

                if (parameter == 'RAJ'):
                    ra = Angle(getattr(self.par_object, 'RAJ'), unit=u.hour)
                    err = getattr(self.par_object, 'RAJerr') / 3600 * factor / n_dim
                    ra_new = ra.deg + err * np.random.uniform(-1., 1.)
                    ra_new = Angle(ra_new, unit=u.deg)
                    setattr(self.par_object, parameter, str(ra_new.to_string(unit=u.hour, sep=':', precision=10)))

                elif (parameter == 'DECJ'):
                    dec = Angle(getattr(self.par_object, 'DECJ'), unit=u.deg)
                    err = getattr(self.par_object, 'DECJerr') / 3600 * factor / n_dim
                    dec_new = dec.deg + err * np.random.uniform(-1., 1.)
                    dec_new = Angle(dec_new, unit=u.deg)
                    setattr(self.par_object, parameter, str(dec_new.to_string(unit=u.deg, sep=':', precision=10)))

                elif (parameter == 'SINI'):
                    cosi = np.sqrt(1 - getattr(self.par_object, 'SINI')**2)
                    err = getattr(self.par_object, 'SINIerr') * factor / n_dim
                    cosi +=  err * np.random.uniform(-1., 1.)
                    setattr(self.par_object, parameter, np.sqrt(1 - cosi**2))

                else:
                    setattr(self.par_object, parameter, value + err * np.random.uniform(-1., 1.))

            setattr(self.par_object, parameter + 'flag', 0)

    def step_size(self, all_sizes=True, parameter_name=None, size=None):

        if all_sizes:

            for parameter in self.parameters:
                self.step_sizes[parameter] = getattr(self.par_object, parameter + 'err')

        else:

            if (parameter_name not in self.parameters):
                sys.exit("Error: option parameter_name not in list of allowed parameters!")

            self.step_sizes[parameter_name] = size
