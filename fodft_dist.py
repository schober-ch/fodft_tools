#!/usr/bin/python
from copy import deepcopy
from fodft_tools import *
import argparse
import os
import sys, traceback
from ase import Atoms

parser = argparse.ArgumentParser(description="Get parameters for fodft")

parser.add_argument('filename', nargs='+', help='Geometry file with the dimer or a list(blob)')#, dest='filename')
parser.add_argument('-e, --extension', help='Format of the geometry file, if not .xyz', dest='fformat', metavar='FORMAT', default='xyz')
parser.add_argument('-d, --dir', help='-d = subfoldername, will create project files there', dest='dir', default='./')
#parser.add_argument('-f, --full', help='Create inputs for basic and polarized fodft', dest='full', action="store_true")
#parser.add_argument('-c, --cubes', help="Automatically adds cube output command for guessed states", dest="cubes", action="store_true")
parser.add_argument('-a, --automagic', help="Tries to find fragments by a clustering algorithm. Check the result carefully! See also:--cutoff", dest="magic", action="store_true")
parser.add_argument('--cutoff', help="Optional: Defines the cutoff for the clustering algorithm. Works best with values larger than inter-molecular distances and smaller than inter-molecular distances! Default is 1.7 for C-C-bonds!", dest="cutoff", type=float) 
#parser.add_argument('-o, --orbitals', help="FO-Orbitals for determination of matrix elements (form: state1 state2 range1 range2)", dest="orbitals", nargs='+', type=int)
parser.add_argument('-i, --image', help="If more than one geometry in .xyz (e.g. a trajectory), which image to choose? Default: Last structure", type=int, dest="image", default=0)
# additional, optinal arguments for the fodft-class
#parser.add_argument('-

arguments = parser.parse_args()

filename = arguments.filename
fformat = arguments.fformat

for file in filename:
    loop = True
    count = 0
    with open(file, "r") as f:
        raw = f.readlines()
    n_atoms = int(raw[0])

    while loop == True:
        try:
            raw[(2+n_atoms)*count]
            count += 1
        except:
            loop = False

    for i in range(count):
        system = fo_aims(file, i, fformat)

        print("Got all information, now create the fragments!")
        if arguments.magic:
            if arguments.cutoff:
                system.magic_cutoff = arguments.cutoff

            system.magic_fragmentation()
        else:
            system.create_fragments()

        print("Now creating the input files!")
        system.write_geom_only()
        print("Done.")
