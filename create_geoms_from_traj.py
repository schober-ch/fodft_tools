#!/usr/bin/python

## takes a xyz with a trajectory (or different distances ..) and writes out a geometry.in for each. HACKZZ!

from ase.io import read,write
import sys

filename = sys.argv[1]

dist = [5.0, 4.5, 4.0, 3.5]

molecules = {}
for num, i in enumerate(dist):
    molecules[i] = read(filename, format="xyz", index=num)


for item in molecules.iteritems():
    write("geometry.in-{0}".format(item[0]), item[1], format="aims")
