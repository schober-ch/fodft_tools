#!/usr/bin/python
import os
import sys
from copy import deepcopy
import ase
from ase.io import read, write
from ase.calculators.aims import Aims

# general stuff

#species_path = "/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/tight"

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

def write_fo(atoms, step):
    # define basic parameters
    if step is "fo":
        str_fo_dft = "final"
    else:
        str_fo_dft = "fragment"
    print("write_fo: {0}, {1}".format(step, str_fo_dft))
    calc = Aims(xc="blyp",
                spin="collinear",
                occupation_type="gaussian 0.01",
                mixer="pulay",
                n_max_pulay="10",
                charge_mix_param="0.5",
                sc_accuracy_rho=1E-4,
                sc_accuracy_eev=1E-2,
                sc_accuracy_etot=1E-6,
                sc_iter_limit=100,
                fo_dft=str_fo_dft,
                default_initial_moment=0,
                relativistic="none",
                species_dir="/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/" + basis_type)

    atoms.set_calculator(calc)

    if str_fo_dft is "fragment":
        atoms.calc.set(fo_options=step)
        if fo_emb is "y":
            atoms.calc.set(fo_embedding=".true.")

    if step is "fo":
        atoms.calc.set(fo_orbitals="x x 1 1 hole")
        atoms.calc.set(packed_matrix_format="none")

    atoms.calc.set(charge=charges[step])
    #atoms.calc.set(default_initial_moment=abs(charges[step])
    

    atoms.calc.write_control(atoms, os.path.join(step, "control.in"))
    atoms.calc.write_species(atoms, os.path.join(step, "control.in"))
                

############################################################################
###################################### CODE ################################
############################################################################

dimer_file = sys.argv[1]

try:
    basis_type = sys.argv[2]
except:
    basis_type = "tight"

print("############################################")
geom_format = raw_input("Format of coord-file (default: xyz): ")

fo_emb = raw_input("Use fo_embedding? [y/n, default: n] ")
if fo_emb == "":
    fo_emb = "n"

charges = {"frag1": 0, "frag2": 0, "fo":0}
t_charges = raw_input("Charges on [frag1], [frag2] (Input-example: +1 0): ")

charges["frag1"] = float(t_charges.strip().split()[0])
charges["frag2"] = float(t_charges.strip().split()[1])
charges["fo"] = charges["frag1"] + charges["frag2"]

if geom_format is "":
    geom_format = "xyz"

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

# create dictionary of atom objects
atom_obj = {"frag1": frag1,
            "frag2": frag2,
            "fo": concat_dimer}

for step in folders:
    write_fo(atom_obj[step], step)


