#! /usr/bin/python

import numpy as np
import sys

class Profile():
    """
    A class that generates a simulated profile and stores its logistical information 
    as attributes.
    """

    def __init__(self, ncomp=1, nbins=1024, fwhm_main=0.05, saveAllComps=False):
        """
        A function that generates the fake (total-intensity_ profile.
        """

        setattr(self, 'n_components', ncomp)
        setattr(self, 'n_bins', nbins)

        pulse_phase = np.linspace(0, 1, num=nbins)
        profile = np.zeros(nbins)

        for component in range(ncomp):
            
            # main pulse always a Gaussian centered at half of pulse phase.
            if (component == 0):
                profile = np.exp(-(pulse_phase - 0.5)**2 / fwhm_main**2)
                if saveAllComps:
                    setattr(self, 'component1', profile)

            # all other components are randomly placed about the main pulse.
            else:
                amplitude = np.random.uniform(0, 0.7)
                width = np.random.uniform(0.01, 0.2)
                location = np.random.uniform(0.4, 0.6)
                comp_profile = amplitude * np.exp(-(pulse_phase - location)**2 / width**2)
                if saveAllComps:
                    setattr(self, 'component' + str(component+1), comp_profile)
                profile += comp_profile

        setattr(self, 'pulse_phase', pulse_phase)
        setattr(self, 'template', profile / np.max(profile))

    def train(self, npulses=1, scintillate=False, noise=False, noise_width=0.2):
        """
        A function that generates a train of pulses, based on the generated template.
        """

        setattr(self, 'n_pulses', npulses)
        pulse_train = np.zeros(self.n_bins * npulses)
        scint_amplitude = 1.

        for num in range(npulses):
        
            if scintillate: 
                scint_amplitude = np.random.uniform(0.0, 1.0)

            pulse_train[num * self.n_bins : (num + 1) * self.n_bins] = scint_amplitude * self.template

        if noise:
            pulse_train += np.random.normal(0, noise_width, size=(npulses * self.n_bins))

        setattr(self, 'pulse_train', pulse_train)
        
