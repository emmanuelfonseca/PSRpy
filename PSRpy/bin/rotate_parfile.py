#! /bin/env python

from PSRpy.parfile import Parfile
import sys

par = Parfile()
par.read(sys.argv[1])
start = par.START["value"]
finish = par.FINISH["value"]
new_epoch = (start + finish) / 2
par.rotate(new_epoch)
par.write()
