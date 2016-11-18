#! /usr/bin/python

import matplotlib.pyplot as plt
import numpy.fft as ft
import numpy as np
import sys

# create train of pulsars.

freq = 0.4
nint = 10.
nsample = 1024 * nint
time = np.linspace(0., nint, num=nsample)
sign = np.sin(2 * np.pi * freq * time)

#plt.plot(time, sign)
#plt.show()

fsig = ft.fft(sign)
ps = fsig.real**2 + fsig.imag**2

timestep = 10. / nsample
ffrq = ft.fftfreq(sign.size, d=timestep)

plt.plot(ffrq, ps)
plt.show()
