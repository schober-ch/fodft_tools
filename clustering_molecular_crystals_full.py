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

parser.add_argument('filename', help='.cif with _chemical_formula_sum!')
parser.add_argument('-d, --cutoff', help='Cutoff-Distance for the clustering algorithm (larger than minimal atomic distance but smaller than molecule-molecule distance', dest='cutoff', default=1.7, type=float)
#parser.add_argument('-n, --number', dest='number_of_molecules', help='Number of molecules in the supercell', type=int)
parser.add_argument('-m, --multiplicator', dest='multiplicator', help='Multiplicator to create supercell from initial cif file. Number depends on the number of unit cells necessary for a full nearest neighbour shell. Script determines largest unit cell vector and tries to adjust the other ratios to generate a cube-like shape. List of 3 ints: 1 1 1. Default: 0', nargs="+", type=int)
parser.add_argument('-c, --central', dest='central_atom', help='Atomic number of a random atom from the central molecule (the one for the neighbour-pairs!', type=int)
parser.add_argument('--com_cutoff', dest='com_cutoff', help='Distance in A for which dimers are generated (based on center of mass, so chain-like molecules will be harder in chain-direction)', default=1000, type=float)

# parse arguments
args = parser.parse_args()

# get data for clustering
supercell = read(args.filename, format="cif")

# only 3 dimensions, therefore hardcoding indices is ok.
box = supercell.get_cell()
box_dim = numpy.array([numpy.linalg.norm(x) for x in box])
box_sort = box_dim.argsort()
vv1 = int(round(box_dim[box_sort[-1]]/box_dim[box_sort[-2]],0))
vv2 = int(round(box_dim[box_sort[-1]]/box_dim[box_sort[-3]],0))

mult_array = [1,1,1]
mult_array[box_sort[-2]] = vv1
mult_array[box_sort[-3]] = vv2

print(mult_array)

full_mult = numpy.array(mult_array) + numpy.array(args.multiplicator)
print(full_mult)

supercell = supercell*full_mult
coordinates = supercell.get_positions()

with open(args.filename, "r") as f:
    cif = f.readlines()
    
for line in cif:
    if "_chemical_formula_sum" in line:
        chem_sum = line.split("_sum")[1].strip().replace(" ", "").strip("'")

clustervector = fclusterdata(coordinates, args.cutoff, criterion="distance")



#print("Number of atoms in each fragment:")
#for i in range(args.number_of_molecules):
    #print(list(clustervector).count(i+1))

frags = []

for i in range(clustervector.max()):
    frags.append(deepcopy(supercell))
    del frags[i][[atom.index for pos, atom in enumerate(frags[i]) if clustervector[pos] != i+1]]

# now get the actual full fragments
final_frags = []
for num, x in enumerate(frags):
    #print([x.get_chemical_formula(), chem_sum])
    if chem_sum in x.get_chemical_formula():
        final_frags.append(frags[num])
    
clean_supercell = final_frags[0]
for i in range(len(final_frags)-1):
    clean_supercell = clean_supercell + final_frags[i+1]

write("clean_supercell.xyz", clean_supercell, format="xyz")    

kdt = KDTree(clean_supercell.get_positions())
nearest = kdt.query(clean_supercell.get_center_of_mass())
pos = clean_supercell[nearest[1]].position

for n, i in enumerate(final_frags):
    for lines in i.get_positions():
        if list(lines) == list(pos):
             central_fragment = n


# now create dimer pairs!
os.mkdir("dimers")
os.chdir("dimers")

print(central_fragment)

# crude way to get nearest neighbour shells easier for chain-like molecules - extrapolate by cubic cell vectors

central_fragment.center(vacuum=1)
cuboid = central_fragment.get_cell()


f = open("center_of_masses.data", "w")

for i in [x for x in range(len(final_frags)) if x != central_fragment]:
    com_dist = numpy.linalg.norm(final_frags[central_fragment].get_center_of_mass() - final_frags[i].get_center_of_mass())
    if float(com_dist) < float(args.com_cutoff):
        write("dimer_{0}_{1}_d-{2}A.xyz".format(central_fragment, i, round(com_dist,2)), (final_frags[central_fragment]+final_frags[i]), format="xyz")
        f.write("COM Central ;{0}; {1}; COM Molec; {2}; {3}; {4}; {5}; Distance; {6}\n".format(central_fragment, final_frags[central_fragment].get_center_of_mass(), i, final_frags[i].get_center_of_mass()[0], final_frags[i].get_center_of_mass()[1], final_frags[i].get_center_of_mass()[2], round(com_dist,2)))
    #f.write("dimer_{0}_{1}_d-{2}A.xyz".format(central_fragment, i, round(com_dist,2)), (final_frags[central_fragment]+final_frags[i]))
    #f.write("COM Fragment {0}: {1}".format(i, final_frags[i].get_center_of_mass()))
    
f.close()
