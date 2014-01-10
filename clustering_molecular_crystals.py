#!/usr/bin/python

from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from copy import deepcopy
import argparse
from ase.io import read, write
import sys, os
import numpy

parser = argparse.ArgumentParser(description="Get parameters to generate neighbour-pairs for molecular crystal xyz")

parser.add_argument('filename', help='XYZ-File (only full molecules in file, no fragments!)')
parser.add_argument('-d, --cutoff', help='Cutoff-Distance for the clustering algorithm (larger than minimal atomic distance but smaller than molecule-molecule distance', dest='cutoff', default=1.7, type=float)
parser.add_argument('-n, --number', dest='number_of_molecules', help='Number of molecules in the supercell', type=int)
parser.add_argument('-c, --central', dest='central_atom', help='Atomic number of a random atom from the central molecule (the one for the neighbour-pairs!', type=int)

args = parser.parse_args()

supercell = read(args.filename, format="xyz")

coordinates = supercell.get_positions()
clustervector = fclusterdata(coordinates, args.cutoff, criterion="distance")

print("The algorithm found {0} molecules!".format(clustervector.max()))
if clustervector.max() == args.number_of_molecules:
    print("Everything is fine.")
else:
    print("Molecules found differ from input! Please check!")
    sys.exit()

print("Number of atoms in each fragment:")
for i in range(args.number_of_molecules):
    print(list(clustervector).count(i+1))

frags = []

for i in range(args.number_of_molecules):
    frags.append(deepcopy(supercell))
    del frags[i][[atom.index for pos, atom in enumerate(frags[i]) if clustervector[pos] != i+1]]

if args.central_atom:
    pos = supercell[args.central_atom].position
else:
    kdt = KDTree(coordinates)
    nearest = kdt.query(supercell.get_center_of_mass())
    pos = supercell[nearest[1]].position

for n, i in enumerate(frags):
    for lines in i.get_positions():
        if list(lines) == list(pos):
             central_fragment = n


# now create dimer pairs!
os.mkdir("dimers")
os.chdir("dimers")

print(central_fragment)
for i in [x for x in range(args.number_of_molecules) if x != central_fragment]:
    com_dist = numpy.linalg.norm(frags[central_fragment].get_center_of_mass() - frags[i].get_center_of_mass())
    #n_dist = cdist(frags[central_fragment].get_positions(), frags[i].get_positions())
    write("dimer_{0}_{1}_d-{2}A.xyz".format(central_fragment, i, round(com_dist,2)), (frags[central_fragment]+frags[i]), format="xyz")
