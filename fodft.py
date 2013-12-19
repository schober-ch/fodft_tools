#!/usr/bin/python

import sys
from fodft_tools_oo import *

filename = sys.argv[1]
try:
    fformat = sys.argv[2]
except:
    fformat = "xyz"

system = fo_aims(filename)

# get params for fodft

xc = raw_input("Which XC functional (Default: blyp):" )
if xc is "":
    system.xc = "blyp"
else:
    system.xc = xc

charge =  raw_input("Charges on [frag1], [frag2] (Input-example: +1 0]): ")
if charge is "":
    system.charge_in = [+1, 0]
else:
    cc = charge.strip("[]").split()
    system.charge_in = cc 

embedding = raw_input("Use embedding? [.true. / .false.] Default no: ")
if embedding is "":
    system.embedding = ".false."
elif embedding is ".true." or "y":
    system.embedding = ".true."

print("Specify basis set, following are possible:")
print(system.avail_species.keys())
basis = raw_input("Default tight. Choice: ")

if basis is "":
    system.species = "tight"
else:
    system.species = basis

#print(system.charge_in)
#print(system.species)
#print(system.embedding)
#print(system.charges)
fotype = raw_input("FO_Type, hole or elec (default: hole)?: ")
if fotype is "":
    system.fo_type = "hole"
else:
    system.fo_type = fotype

print("Got all information, now create the fragments!")
system.create_fragments()

print("Now creating the input files!")
system.write_inputs()

print("Done.")
