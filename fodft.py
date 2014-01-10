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
parser.add_argument('-f, --full', help='Create inputs for basic and polarized fodft', dest='full', action="store_true")
parser.add_argument('-c, --cubes', help="Automatically adds cube output command for guessed states", dest="cubes", action="store_true")
parser.add_argument('-a, --automagic', help="Tries to find fragments by a clustering algorithm. Check the result carefully! See also:--cutoff", dest="magic", action="store_true")
parser.add_argument('--cutoff', help="Optional: Defines the cutoff for the clustering algorithm. Works best with values larger than inter-molecular distances and smaller than inter-molecular distances! Default is 1.7 for C-C-bonds!", dest="cutoff", type=float) 
parser.add_argument('-o, --orbitals', help="FO-Orbitals for determination of matrix elements (form: state1 state2 range1 range2)", dest="orbitals", nargs='+', type=int)

# additional, optinal arguments for the fodft-class
#parser.add_argument('-

arguments = parser.parse_args()

filename = arguments.filename
fformat = arguments.fformat

#if len(filename) > 1 and arguments.dir:
    #print("Using one folder and more than one input doesn't work! Bye!")
    #sys.exit()
example = fo_aims(Atoms('CO', positions=[(0, 0, 0), (0, 0, 1)]))
arg_dict = {
            "xc" :        ["Which XC functional (Default: blyp): ", "blyp"],
            "charge_in" :    ["Charges on [frag1], [frag2] (Default: +1 0]): ", "+1 0"],
            "embedding" : ["Use embedding? [y/n] (Default: no): ", ".false."],
            "species" :     ["Specify basis set, available options: \n\n {0} \n\n(Default: tight). Please choose: ".format(example.avail_species.keys()), "tight"],
            "fo_type" :    ["FO_Type, hole or elec (Default: hole): ", "hole"],
            }

# get params for fodft
arg_dict_values = deepcopy(arg_dict)

# First, get user input
for item in arg_dict:
    arg_dict_values[item][1] = raw_input("{0}".format(arg_dict[item][0]))

# Fill up with defaults were no user input exists
for item in arg_dict_values:
    if arg_dict_values[item][1] is "":
        arg_dict_values[item][1] = arg_dict[item][1]

# Special post processing of inputs
if arg_dict_values['embedding'][1] == "y":
    arg_dict_values['embedding'][1] = ".true."
elif arg_dict_values['embedding'][1] == "n":
    arg_dict_values['embedding'][1] = ".false."

for file in filename:
    system = fo_aims(file, fformat)

    if len(filename) > 1:
        dirname = file.rpartition(".")[0]
        if dirname == filename:
            dirname = "fo" + dirname
        arguments.dir = dirname

    if arguments.dir == "./":
        print("Creating files in current working directory ({0})".format(os.getcwd()))
    else:
        try:
            cwd_start = os.getcwd()
            os.mkdir(arguments.dir)
            os.chdir(arguments.dir)
            print("Creating files in {0}!".format(arguments.dir))
        except:
            print("Error when creating folder {0}:".format(arguments.dir))
            raise

    print("Creating basic and embedded input: {0}".format(arguments.full))

    # Set the values. 
    system.aims_params['xc'] = arg_dict_values['xc'][1]
    system.charge_in = arg_dict_values['charge_in'][1].strip("[]").split()
    system.embedding = arg_dict_values['embedding'][1]
    system.species = arg_dict_values['species'][1]
    system.fo_type = arg_dict_values['fo_type'][1]

    print("Got all information, now create the fragments!")
    if arguments.magic:
        if arguments.cutoff:
            system.magic_cutoff = arguments.cutoff

        system.magic_fragmentation()
    else:
        system.create_fragments()

    if arguments.orbitals:
        system.frontiers = arguments.orbitals[0:2]
        system.fo_range = arguments.orbitals[2:4]
        system.update_calculators()
    
    if arguments.cubes is True:
        system.set_cube_files()

    if arguments.full is True:
        print("Now creating input files for basic fo_dft...")

        os.mkdir("basic")
        os.mkdir("embedded")

        cwd = os.getcwd()
        os.chdir("basic")
        #print("Now creating the input files!")
        system.write_inputs()
        os.chdir(cwd)
        os.chdir("embedded")
        print("Now creating input files for embedded fo_dft...")
        system.embedding = ".true."
        system.update_calculators()
        system.set_cube_files()
        system.write_inputs()
    else:
        print("Now creating the input files!")
        system.write_inputs()
    print("Done.")
    os.chdir(cwd_start)
