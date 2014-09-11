#!/sw/bin/python -O

from numpy import *

class Profile(object):
    def __init__(self, filename, suffix="",reticon=False):
        f = open(filename+suffix)
        list = f.readline().split()
        self.frequency = float(list[5])
        self.bins = int(list[7])
        self.name = list[10]
        self.bin = zeros(self.bins, dtype=float)
        self.power = zeros(self.bins, dtype=float)
        self.q = zeros(self.bins, dtype=float)
        self.u = zeros(self.bins, dtype=float)
        self.linear = zeros(self.bins, dtype=float)
        self.circular = zeros(self.bins, dtype=float)
        self.angle = zeros(self.bins, dtype=float)
        self.angle_error = zeros(self.bins, dtype=float)
        for i in xrange(self.bins):
            list = f.readline().split()
            self.bin[i] = int(list[0])-1
            self.power[i] = float(list[1])
            self.q[i] = float(list[2])
            self.u[i] = float(list[3])
            self.circular[i] = float(list[4])
            self.linear[i] = float(list[5])
            self.angle[i] = float(list[6])
            self.angle_error[i] = float(list[7])
	if reticon:
            self.linear = self.u
            self.circular = self.q
