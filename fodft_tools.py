#!/usr/bin/python
#

import os
import sys
from copy import deepcopy
import ase
from ase.io import read,write
from ase.calculators.aims import Aims
from ase.atoms import atomic_numbers
from ase.calculators.cpmd import CPMD

# global parameters
spec_path = "/data/schober/code/fhiaims_develop/fhiaims_supporting_work/species_defaults/"
avail_species = {"light" : "light",
                 "tight" : "tight",
                 "cc.3"  : "non-standard/NAO-VCC-nZ/NAO-VCC-3Z",
                 "cc.4"  : "non-standard/NAO-VCC-nZ/NAO-VCC-4Z",
                 "cc.4"  : "non-standard/NAO-VCC-nZ/NAO-VCC-5Z",
                 "tight.ext" : "tight.ext",
                 "cc.3.ext"  : "non-standard/NAO-VCC-nZ/NAO-VCC-3Z.ext"
                 }

class fodft:
    """ Functionality for fragment orbital dft calculations. This class collects all the methods necessary for FODFT calculations, regardless of the used QM program. """

    def __init__(self, dimer, fformat="xyz"):
        # Atom objects
        self.dimer = None
        self.frag1 = None
        self.frag2 = None

        # General FO-DFT parameters
        self.charge_in = [0, 0]
        self.charges = {"frag1": 0, "frag2": 0, "fo":0}
        self.initial_moments = None
        self.fo_type = "hole"
        self.frontiers = [1, 1]
        self.fo_range = [2, 2]
        try:
            dimer.get_name()
            self.dimer = dimer
        except:
            self.read_dimer(dimer, fformat)

    def read_dimer(self, filename, fformat):
        """ Read the dimer file with a given file format and call create_fragments after that. """
        
        try:
            self.dimer = read(filename, format=fformat)
        except:
            print("ERROR: The file {0} was not found or the format {1} not known!".format(filename, fformat))

    def create_fragments(self):
        """ This function creates fragments from a dimer file by letting you visually select the fragments. Format: [charge1, charge2]."""

        self.frag1 = deepcopy(self.dimer)
        self.frag2 = deepcopy(self.dimer)

        print("Edit fragment 1")
        self.frag1.edit()
        print("Edit fragment 2")
        self.frag2.edit()
        print("Fragments created.")

        check = (self.frag1 + self.frag2).get_chemical_formula() == self.dimer.get_chemical_formula()
        if check is True:
            print("Created dimers, checked consistency.")
        else:
            print("ERROR: Number of atoms of frag1 + frag2 is not the same as dimer! Please create fragments again!") 
            self.create_fragments()

        self.__set_charges__()
        self.__get_frontiers__()

    def __set_charges__(self):
        """ Sets the charges for fragments and dimer according to the input.
    Usage: 
        obj.set_charges([+1, 0])
    Example:
        Fragment 1 has   +1
        Fragment 2 has    0
        Dimer must have  +1"""

        # Initialize standard values
        self.charges["frag1"] = float(self.charge_in[0])
        self.charges["frag2"] = float(self.charge_in[1])
        self.charges["fo"] = self.charges["frag1"] + self.charges["frag2"]

        # update initial moments after charge changed!
        self.__get_initial_moments__()

    def __get_frontiers__(self):
        """ Extracts the HOMOs for both fragments in their _neutral_ state to be used for determination of the orbitals of interest.""" 
        self.frontiers[0] = self.frag1.get_atomic_numbers().sum()
        self.frontiers[1] = self.frag2.get_atomic_numbers().sum()

    def __get_initial_moments__(self):
        """ Calculates initial moments from number of electrons and charge. Only works for low spin systems! ALWAYS check before submitting! """
        total_f1 = self.frag1.get_atomic_numbers().sum() - self.charges["frag1"]
        total_f2 = self.frag2.get_atomic_numbers().sum() - self.charges["frag2"]
        total_fo = self.dimer.get_atomic_numbers().sum() - self.charges["fo"]
        
        # now, div by 2 or not? (I know, bad way..)
        self.initial_moments = [total_f1%2, total_f2%2, total_fo%2]

class fo_aims(fodft):
    """ Class for fragment orbital dft calculations with FHI_aims. """
        
    def __init__(self, dimer, fformat="xyz"):
        fodft.__init__(self, dimer, fformat)
        
        self.avail_species = avail_species 
        self.species = "tight"
        self.embedding = ".false."

        # Aims standard-params
        self.aims_params = {
                    "xc"     : "blyp",
                    "spin"  : "collinear",
                    "occupation_type"   : "gaussian 0.01",
                    "mixer"     : "pulay",
                    "n_max_pulay"   : "10",
                    "charge_mix_param" : "0.5",
                    "sc_accuracy_rho"   : "1E-4",
                    "sc_accuracy_eev"   : "1E-2",
                    "sc_accuracy_etot"  : "1E-6",
                    "relativistic"  : "none",
                    "species_dir" : os.path.join(spec_path, avail_species[self.species])
                    }
           # add a calculator
        #self.update_calculators()
    def create_fragments(self):
        fodft.create_fragments(self)
        self.update_calculators()

    def update_calculators(self):
        self.__set_charges__()

        self.aims_params['species_dir'] = os.path.join(spec_path, avail_species[self.species])
        self.dimer.set_calculator(Aims(**self.aims_params))
        self.frag1.set_calculator(Aims(**self.aims_params))
        self.frag2.set_calculator(Aims(**self.aims_params))

        # params for fragments:
        self.frag1.calc.set(default_initial_moment=self.initial_moments[0], 
                            fo_dft="fragment",
                            fo_options="frag1",
                            fo_embedding=self.embedding,
                            charge=self.charges['frag1'])
        self.frag2.calc.set(default_initial_moment=self.initial_moments[1], 
                            fo_dft="fragment", 
                            fo_options="frag2",
                            fo_embedding=self.embedding,
                            charge=self.charges['frag2'])

        self.dimer.calc.set(default_initial_moment=self.initial_moments[2], 
                         fo_dft="final", 
                         charge=self.charges['fo'],
                         fo_orbitals="{0} {1} {2} {3} {4}".format(self.frontiers[0]/2, self.frontiers[1]/2, self.fo_range[0], self.fo_range[1], self.fo_type),
                         packed_matrix_format="none")

    def write_inputs(self):

        #fragments
        os.mkdir("frag1")
        os.mkdir("frag2")

        self.frag1.calc.write_control(self.frag1, "frag1/control.in")
        self.frag1.calc.write_species(self.frag1, "frag1/control.in")
        write("frag1/geometry.in", self.frag1, format="aims")
        
        self.frag2.calc.write_control(self.frag2, "frag2/control.in")
        self.frag2.calc.write_species(self.frag2, "frag2/control.in")

        if self.embedding == ".false.":
            write("frag2/geometry.in", self.frag2, format="aims")
        else:
            self.__write_aims_empty__("frag1/geometry.in", self.frag2)     
            self.__write_species_empty__(self.frag2, "frag1/control.in")
            
            write("frag2/tmp_geometry.in", self.frag2, format="aims")
            self.__write_aims_empty__("frag2/geometry.in", self.frag1)
            self.__write_species_empty__(self.frag1, "frag2/control.in")
            with open('frag2/geometry.in', 'a') as outfile:
                with open("frag2/tmp_geometry.in") as infile:
                        outfile.write(infile.read())
            os.remove("frag2/tmp_geometry.in")

        # final step
        self.dimer.calc.write_control(self.dimer, "control.in")
        self.dimer.calc.write_species(self.dimer, "control.in")
        write("geometry.in", self.frag1+self.frag2, format="aims")

    # this is the write_aims from ASE
    def __write_aims_empty__(self, filename, atoms):
        """Method to write FHI-aims geometry files.

        Writes the atoms positions and constraints (only FixAtoms is
        supported at the moment).

        This method was changed from original ASE io/aims.py to allow
        empty sites. 
        """

        if isinstance(atoms, (list, tuple)):
            if len(atoms) > 1:
                raise RuntimeError("Don't know how to save more than "+
                                   "one image to FHI-aims input")
            else:
                atoms = atoms[0]

        fd = open(filename, 'a')
        fd.write('#=======================================================\n')
        fd.write('#FHI-aims file: '+filename+'\n')
        fd.write('#Created using the Atomic Simulation Environment (ASE)\n')
        fd.write('#=======================================================\n')
        for i, atom in enumerate(atoms):
            fd.write('empty ')
            for pos in atom.position:
                fd.write('%16.16f ' % pos)
            fd.write(atom.symbol+"_e")
            fd.write('\n')

    # Also modified from ASE aims calculator

    def __write_species_empty__(self, atoms, filename):
        species_path = atoms.calc.parameters.get('species_dir')

        control = open(filename, 'a')
        symbols = atoms.get_chemical_symbols()
        symbols2 = []
        for n, symbol in enumerate(symbols):
            if symbol not in symbols2:
                symbols2.append(symbol)
        for symbol in symbols2:
            fd = os.path.join(species_path, '%02i_%s_default' %
                              (atomic_numbers[symbol], symbol))
            for line in open(fd, 'r'):
                if 'species' in line:
                    control.write(line.strip()+"_e")
                    control.write("\n")
                elif '#  "First tier"' in line or '# (sp) correlation set' in line:
                    control.write("include_min_basis .false.")
                    control.write("\n")
                    break
                else:
                    control.write(line)
                    
        control.close()

class fo_cpmd(fodft):
    """ Wrapper for old cpmd io/calculator. Need to integrate into this later on..."""
   
    def __init__(self, dimer, fformat="xyz", charges=1):
        fodft.__init__(self, dimer, fformat="xyz")

        self.dimer.calc = CPMD(self.dimer)
        self.dimer.calc.center(vacuum=4)
        self.dimer.calc.set(charge=charges)

    def create_fragments(self, fo_files=True):
        self.dimer.calc.create_fragments(self.dimer, fo_files=True)


