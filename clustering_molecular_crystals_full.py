#!/usr/bin/python

from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from copy import deepcopy
import argparse
from ase.io import read, write
from ase.visualize import view
import sys, os
import numpy as np
import scipy
from sklearn.cluster import DBSCAN

parser = argparse.ArgumentParser(description="Get parameters to generate neighbour-pairs for molecular crystal xyz")

parser.add_argument('filename', help='.cif with _chemical_formula_sum!')
parser.add_argument('-d, --cutoff', help='Cutoff-Distance for the clustering algorithm (larger than minimal atomic distance but smaller than molecule-molecule distance', dest='cutoff', default=1.7, type=float)
#parser.add_argument('-n, --number', dest='number_of_molecules', help='Number of molecules in the supercell', type=int)
parser.add_argument('-m, --multiplicator', dest='multiplicator', help='Multiplicator to create supercell from initial cif file. Number depends on the number of unit cells necessary for a full nearest neighbour shell. Script determines largest unit cell vector and tries to adjust the other ratios to generate a cube-like shape. Default: 2', type=int, default=2)
parser.add_argument('-c, --central', dest='central_atom', help='Atomic number of a random atom from the central molecule (the one for the neighbour-pairs!', type=int)
parser.add_argument('--com_cutoff', dest='com_cutoff', help='Distance in A for which dimers are generated (based on center of mass, so chain-like molecules will be harder in chain-direction)', default=1000, type=float)
parser.add_argument('--dry', dest='dryrun', help='If choose, only the steps prior to the fragmentation are done and the resulting supercell is shown.', action="store_true")

# parse arguments
args = parser.parse_args()

# get data for clustering
supercell = read(args.filename, format="cif")

# only 3 dimensions, therefore hardcoding indices is ok.
box = supercell.get_cell()
box_dim = np.array([np.linalg.norm(x) for x in box])
box_sort = box_dim.argsort()
vv1 = int(round(box_dim[box_sort[-1]]/box_dim[box_sort[-2]],0))
vv2 = int(round(box_dim[box_sort[-1]]/box_dim[box_sort[-3]],0))

mult_array = [1,1,1]
mult_array[box_sort[-2]] = vv1
mult_array[box_sort[-3]] = vv2

full_mult = np.array(mult_array) * args.multiplicator
print("Unit cell is multiplied by {0}.".format(full_mult))
supercell = supercell*full_mult

if args.dryrun is True:
    view(supercell)
    sys.exit()

coordinates = supercell.get_positions()

with open(args.filename, "r") as f:
    cif = f.readlines()

try:
    for line in cif:
        if "_chemical_formula_sum" in line:
            chem_sum = line.split("_sum")[1].strip().replace(" ", "").strip("'")
except:
    print("No valid chemical formula was found in the cif file. Please check!")
    raise

os.mkdir("dimers")
os.chdir("dimers")

print("Clustering step started...")
#clustervector = fclusterdata(coordinates, args.cutoff, criterion="distance")
db = DBSCAN(eps=args.cutoff, min_samples=1).fit(coordinates)
clustervector = db.labels_.astype(int)
print("Clustering step finished...")

molecule = []
counts = scipy.stats.itemfreq(clustervector)
number_of_atoms = counts[:,1].max()
for item in counts:
    if item[1] == number_of_atoms:
        molecule.append(int(item[0]))

for num, item in enumerate(clustervector):
    if item not in molecule:
        clustervector[num] = molecule[0]

counts = scipy.stats.itemfreq(clustervector)
print("Preparing fragment objects...")

frags = []
for i in range(len(counts[:,0])):
    frags.append(deepcopy(supercell))

for num, i in enumerate(frags):
    del frags[num][[atom.index for pos, atom in enumerate(frags[num]) if clustervector[pos] != counts[num][0]]]

# now get the actual full fragments
# skipped, not neccessary with counts[:,1].max(). The idea is simple: The largest fragment found should reasonably be the molecule, if not, something went wrong...
##final_frags = deepcopy(frags)
final_frags = []
for num, x in enumerate(frags):
    if len(x.get_positions()) == number_of_atoms:
        final_frags.append(frags[num])
    else:
        write("removed_atoms.xyz", x, format="xyz")

print("Found {0} complete fragments.".format(len(counts[:,0])-1))
print("Sanity check: Used {0} of {1} atoms.".format(sum([len(y) for y in [x.get_positions() for x in frags]]), len(supercell.get_positions())))
print("Done.")


# get cleaned supercell
clean_supercell = final_frags[0]
for i in range(len(final_frags)-1):
    clean_supercell = clean_supercell + final_frags[i+1]

# find the molecule closest to the center of the cell
print("Look for center molecule...")
kdt = KDTree(clean_supercell.get_positions())
nearest = kdt.query(clean_supercell.get_center_of_mass())
pos = clean_supercell[nearest[1]].position

for n, i in enumerate(final_frags):
    for lines in i.get_positions():
        if list(lines) == list(pos):
             central_fragment = n

print("Found center molecule, {0}".format(central_fragment))

# now create dimer pairs!

neighbour_frags = []

# use copy of final frags here to avoid handling the central fragment twice
search_frags = deepcopy(final_frags)
search_frags.pop(central_fragment)

print("Started search for nearest neighbours...")
for fragment in search_frags:
    dist = scipy.spatial.distance.cdist(final_frags[central_fragment].get_positions(), fragment.get_positions())
    
    neighbours = {}
    for num, pairs in enumerate(dist):
        if pairs.min() < float(args.com_cutoff):
            neighbours[num] = True
        else:
            neighbours[num] = False

    if any(x is True for x in neighbours.itervalues()):
        neighbour_frags.append(fragment)

neighbour_supercell = final_frags[central_fragment]
for i in range(len(neighbour_frags)):
    neighbour_supercell = neighbour_supercell + neighbour_frags[i]

print("Finished search...")

# write the clean supercell and the neighbour cell to files
write("neighbour_supercell.xyz", neighbour_supercell, format="xyz")
write("clean_supercell.xyz", clean_supercell, format="xyz")    

f = open("center_of_masses.data", "w")

for i in range(len(neighbour_frags)):
    com_dist = np.linalg.norm(final_frags[central_fragment].get_center_of_mass() - neighbour_frags[i].get_center_of_mass())
    write("dimer_{0}_{1}_d-{2}A.xyz".format(central_fragment, i, round(com_dist,2)), (final_frags[central_fragment]+neighbour_frags[i]), format="xyz")
    f.write("COM Central ;{0}; {1}; COM Molec; {2}; {3}; {4}; {5}; Distance; {6}\n".format(central_fragment, final_frags[central_fragment].get_center_of_mass(), i, neighbour_frags[i].get_center_of_mass()[0], neighbour_frags[i].get_center_of_mass()[1], neighbour_frags[i].get_center_of_mass()[2], round(com_dist,2)))
    
f.close()
