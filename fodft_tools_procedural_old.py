#!/usr/bin/python
import os
import sys
from copy import deepcopy
import ase
from ase.io import read, write
from ase.calculators.aims import Aims

# general stuff

def create_fragments(atoms, fo_files=None):
    frag1 = deepcopy(atoms)
    frag2 = deepcopy(atoms)

    print("Edit fragment 1")
    frag1.edit()

    print("Edit fragment 2")
    frag2.edit()

    print("Fragments created.")

    return frag1, frag2

def write_fo(atoms, step, orbitals):
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
                #fo_dft=str_fo_dft,
                default_initial_moment=abs(charges[step]),
                relativistic="none",
                species_dir="/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/" + basis_type)

    atoms.set_calculator(calc)

    if basis_type.split(".")[0] == "cc":
        #if basis_type.split(".")[1] is "3" or "4" or "5":
        atoms.calc.set(species_dir="/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/non-standard/NAO-VCC-nZ/NAO-VCC-{0}Z".format(basis_type.split(".")[1]))
        #else:
            #print("This cc-basis, {0}, does not exist!".format(basis_type))

    if str_fo_dft is "fragment":
        atoms.calc.set(fo_options=step)
        if fo_emb is "y":
            #atoms.calc.set(fo_embedding=".true.")

    if step is "fo":
        
        atoms.calc.set(fo_orbitals="{0} {1} 3 3 hole".format(orbitals[0]/2, orbitals[1]/2))
        atoms.calc.set(packed_matrix_format="none")
    atoms.calc.set(charge=charges[step])
    #atoms.calc.set(default_initial_moment=abs(charges[step])
    

    if step is not "fo":
        atoms.calc.write_control(atoms, os.path.join(step, "control.in"))
        atoms.calc.write_species(atoms, os.path.join(step, "control.in"))
        if fo_emb is "y":
            # need to change manually!
            atoms.calc.write_species(atoms, os.path.join(step, "control.in"))

    elif step is "fo":
        atoms.calc.write_control(atoms, "control.in")
        atoms.calc.write_species(atoms, "control.in")

############################################################################
###################################### CODE ################################
############################################################################

dimer_file = sys.argv[1]

# debug:
#dimer_file = "two_thiophene.xyz"
try:
    basis_type = sys.argv[2]
#except:
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

# get electron number for fo dft (HOMO without charges)
orbitals = []
orbitals.append(frag1.get_atomic_numbers().sum())
orbitals.append(frag2.get_atomic_numbers().sum())

print("""Orbital numbers (HOMO) without charges for frag1: {0} 
... and frag2: {1}""".format(orbitals[0]/2, orbitals[1]/2))

#create frag directories
folders = ["frag1", "frag2", "fo"]

for folder in folders:
    os.mkdir(folder)

if fo_emb == "y":
    write("frag1/geometry.in", frag1, format="aims")
    write("frag1/geometry.in", frag2, format="aims")
    write("frag2/geometry.in", frag1, format="aims")
    write("frag2/geometry.in", frag2, format="aims")
else:
    write("frag1/geometry.in", frag1, format="aims")
    write("frag2/geometry.in", frag2, format="aims")

concat_dimer = frag1 + frag2

write("geometry.in", concat_dimer, format="aims")

# create dictionary of atom objects
atom_obj = {"frag1": frag1,
            "frag2": frag2,
            "fo": concat_dimer}

for step in folders:
    write_fo(atom_obj[step], step, orbitals)


