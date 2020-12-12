#! /usr/bin/python

import numpy as np
import sys

class MCMC():

    def __init__(self, xdata, ydata, uncertainties, parameter_names=[]):
        """
        Initialize sampler class.
        """

        self.xdata = xdata
        self.ydata = ydata
        self.n_params = len(parameter_names)
        self.parameters = parameter_names
        self.uncertainties = uncertainties
        self.priors = {}
        self.posteriors = {}
        self.step_sizes = {}

    def log_likelihood(self, parameter_values):
        """
        Computes the likelihood.

        Inputs
        ------

        model_function : callable
            Function that computes model from input parameters.
        """

        model_ydata = self.model(self.xdata, p=parameter_values)
        chisq = np.sum((model_ydata - self.ydata)**2 / self.uncertainties)
        return -0.5 * chisq

    def model(self, model_function):
        """
        Computes chi-squared goodness of fit statistic.

        Inputs
        ------

        model_function : callable
            Function that computes model from input parameters.
        """

        self.model = model_function

    def plot(self, n_bins=100):
        """
        Plots MCMC posterior data.
        """

        from matplotlib.font_manager import FontProperties
        import matplotlib.pyplot as plt

        font = FontProperties()
        font.set_name('serif')
        plt.style.use('classic')

        # plot chains versus iterations.
        plt.figure(1)

        for num, parameter in zip(range(self.n_params), self.parameters):
            plt.subplot(self.n_params, 1, num + 1)
            plt.plot(self.posteriors[parameter], 'b-')

        plt.figure(2)

        from matplotlib.gridspec import GridSpec
        gs = GridSpec(self.n_params, self.n_params)

        count = 0

        for num1, parameter1 in zip(range(self.n_params), self.parameters):
            for num2, parameter2 in zip(range(self.n_params), self.parameters):
                if (num1 == num2):
                    plt.subplot(gs[num1, num2])
                    plt.hist(self.posteriors[parameter1], bins=n_bins, histtype='step')
                    plt.xticks([])
                    plt.yticks([])
                elif (num1 > num2):
                    plt.subplot(gs[num1, num2])
                    hist2D, xedges, yedges = np.histogram2d(self.posteriors[parameter2], self.posteriors[parameter1], bins=50)
                    dx = xedges[1] - xedges[0]
                    dy = yedges[1] - yedges[0]
                    x = np.linspace(xedges[0] + dx, xedges[len(xedges)-2] + dx, num=len(xedges)-1)
                    y = np.linspace(yedges[0] + dy, yedges[len(yedges)-2] + dy, num=len(yedges)-1)
                    plt.pcolormesh(x, y, hist2D, cmap='Blues')
                else:
                    pass
                count += 1

        plt.show()

    def prior(self, prior_type, all_priors=True, parameter_name=None):
        """
        Set type of prior distribution for one or all parameters.

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

    def run(self, parameter_values, burn=500, n_iterations=1000, thin=0):
        """
        Runs the MCMC sampler.

        Inputs
        ------

        parameter_values : array_like
            An array of values for the parameters to be sampled.

        n_iterations : int
            Number of iterations for sampler. Default is 1000.

        thin : int
            Thinning factor. Default is 0.
        """

        samples = np.zeros((n_iterations, self.n_params))
        success_count = 0
        current_par = parameter_values

        for iternum in range(n_iterations):

            proposal_par = self.step(current_par)
            ll_current = self.log_likelihood(current_par)
            ll_proposal = self.log_likelihood(proposal_par)
            log_prior = np.log(np.random.uniform(0., 1.))

            if (log_prior < np.min([0., ll_proposal - ll_current])):
                current_par = proposal_par
                success_count += 1
                samples[iternum, :] = proposal_par

            else:
                samples[iternum, :] = current_par

        self.accepted = success_count

        for num, parameter in zip(range(self.n_params), self.parameters):
            self.posteriors[parameter] = samples[burn:, num]

    def step(self, parameter_values):
        """
        Steps a set of parameters, given information on prior.
        """

        new_parameters = np.zeros(self.n_params)        
        count = 0

        for parameter, val in zip(self.parameters, parameter_values):
            if (self.priors[parameter] == 'uniform'):
                new_parameters[count] = val + np.random.uniform(-1, 1) * self.step_sizes[parameter]
            elif (self.priors[parameter] == 'normal'):
                new_parameters[count] = val + np.random.normal(0, 1) * self.step_sizes[parameter]
            count += 1

        return new_parameters

    def step_size(self, size, all_sizes=True, parameter_name=None):

        if all_sizes:

            for parameter in self.parameters:
                self.step_sizes[parameter] = size

        else:

            if (parameter_name not in self.parameters):
                sys.exit("Error: option parameter_name not in list of allowed parameters!")

            self.step_sizes[parameter_name] = size
