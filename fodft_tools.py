import os
import sys
from copy import deepcopy
import ase
from ase.io import read, write
from ase.calculators.aims import Aims

# general stuff

species_path = "/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/tight"

def create_fragments(atoms, fo_files=None):
    """ Function to create fragments for FO-DFT from complete cell
        
        Usage: 
        
            frag1, frag2 = atoms.calc.create_fragments(atoms)
        or                        
            atoms.calc.create_fragments(atoms, fo_files=True)

        Latter will invoke atoms.calc.fodft and create the inputs
               
               """

    frag1 = deepcopy(atoms)
    frag2 = deepcopy(atoms)

    print("Edit fragment 1")
    frag1.edit()

    print("Edit fragment 2")
    frag2.edit()

    print("Fragments created.")

    return frag1, frag2



dimer_file = sys.argv[1]

print("############################################")
geom_format = raw_input("Format of coord-file (default: xyz): ")

dimer = read(dimer_file, format=geom_format)

#aimscalc = Aims(xc='PBE')

#dimer.set_calculator(aimscalc)

frag1, frag2 = create_fragments(dimer)

currdir = os.getcwd()

#create frag directories
folders = ["frag1", "frag2", "fo"]

for folder in folders:
    os.mkdir(folder)

write("frag1/geometry.in", frag1, format="aims")
write("frag2/geometry.in", frag2, format="aims")

concat_dimer = frag1 + frag2

write("fo/geometry.in", concat_dimer, format="aims")

for step in folders:
    write_fo(step)


