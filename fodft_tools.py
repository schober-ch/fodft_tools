#!/usr/bin/python
#

import os
import sys
from copy import deepcopy
import ase
from ase.io import read,write
from ase.calculators.aims import Aims
from ase.atoms import atomic_numbers

#from scipy.cluster.vq import kmeans,vq # used for the automagic clustering
from scipy.cluster.hierarchy import fclusterdata # used for the automagic clustering, better algorithm than scipy.cluster.vq

#from ase.calculators.cpmd import CPMD

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

    def __init__(self, dimer, fformat="xyz", w_image):
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
        self.magic_cutoff = 1.7 # default for c-c bond distance
        try:
            dimer.get_name()
            self.dimer = dimer
        except:
            self.read_dimer(dimer, fformat, w_image)

    def read_dimer(self, filename, fformat, w_image):
        """ Read the dimer file with a given file format and call create_fragments after that. """
        
        try:
            self.dimer = read(filename, format=fformat, index=w_image)
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

        self.__check_fragments__()

        self.__set_charges__()
        self.__get_frontiers__()

    def magic_fragmentation(self):
        """ This function takes the atom objects and tries to separate two fragments by a k-means-clustering algorithm. Always check the result before relying on those fragmentations!"""

        #hardcoded number of fragments, for now always 2!
        nr_frags = 2
        coordinates = self.dimer.get_positions()
        # 
        #centroids,_ = kmeans(coordinates, nr_frags)
        # assign indices to clusters (bitmask!)
        cluster_indices = fclusterdata(coordinates, self.magic_cutoff, criterion="distance")
        # compress the whole coordinates to fragments 
        #coords_frag1 = np.array(list(itertools.compress(coordinates.tolist(), cluster_indices)))
        # invert the bitmask
        #cluster_indices = cluster_indices ^ 1
        #coords_frag2 = np.array(list(itertools.compress(coordinates.tolist(), cluster_indices)))

        self.frag1 = deepcopy(self.dimer)
        self.frag2 = deepcopy(self.dimer)

        # Now delete the atoms of the other fragment from the object with mighty pythonic list comprehensions!
        del self.frag1[[atom.index for pos, atom in enumerate(self.frag1) if cluster_indices[pos] != 1]]
        del self.frag2[[atom.index for pos, atom in enumerate(self.frag2) if cluster_indices[pos] != 2]]

        print("Finished automatic fragmentation, please remember to check the result!")
        self.__check_fragments__() 

        self.__set_charges__()
        self.__get_frontiers__()

    def __check_fragments__(self):
        """ This function contains all checks for the fragmentation process, such as comparing the total number of atoms. """
        print("Chemical formula for fragment1 and fragment2: {0} | {1}".format(self.frag1.get_chemical_formula(), self.frag2.get_chemical_formula()))

        check = (self.frag1 + self.frag2).get_chemical_formula() == self.dimer.get_chemical_formula()
        if check is not True:
            print("ERROR: Number of atoms of frag1 + frag2 is not the same as dimer! Please create fragments again!")
            sys.exit()
        else:
            print("Created dimers, did a _basic_ check for consistency.") 

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
        self.frontiers[0] = self.frag1.get_atomic_numbers().sum()/2
        self.frontiers[1] = self.frag2.get_atomic_numbers().sum()/2

    def __get_initial_moments__(self):
        """ Calculates initial moments from number of electrons and charge. Only works for low spin systems! ALWAYS check before submitting! """
        total_f1 = self.frag1.get_atomic_numbers().sum() - self.charges["frag1"]
        total_f2 = self.frag2.get_atomic_numbers().sum() - self.charges["frag2"]
        total_fo = self.dimer.get_atomic_numbers().sum() - self.charges["fo"]
        
        # now, div by 2 or not? (I know, bad way..)
        self.initial_moments = [total_f1%2, total_f2%2, total_fo%2]

class fo_aims(fodft):
    """ Class for fragment orbital dft calculations with FHI_aims. """
        
    def __init__(self, dimer, fformat="xyz", w_image):
        fodft.__init__(self, dimer, fformat, w_image)
        
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

    def magic_fragmentation(self):
        fodft.magic_fragmentation(self)
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
                         fo_orbitals="{0} {1} {2} {3} {4}".format(self.frontiers[0], self.frontiers[1], self.fo_range[0], self.fo_range[1], self.fo_type),
                         packed_matrix_format="none")

    def set_cube_files(self):
        """ Method to create cube file input for each fragment """
        
        # its done for both fragments each, might be able to do this nicer..    
        states = range(self.frontiers[0], self.frontiers[0] + self.fo_range[0])
        plots = ["eigenstate {0}".format(x) for x in states]
        cube1 = AimsCube(plots)
        self.frag1.calc.cubes = cube1

        states = range(self.frontiers[1], self.frontiers[1] + self.fo_range[1])
        plots = ["eigenstate {0}".format(x) for x in states]
        cube2 = AimsCube(plots)
        self.frag2.calc.cubes = cube2

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

# also copied from ASE aims calculator, modified for fodft.. 
class AimsCube:
    "Object to ensure the output of cube files, can be attached to Aims object"
    def __init__(self, plots=None):
        """parameters:
        origin, edges, points = same as in the FHI-aims output
        plots: what to print, same names as in FHI-aims """

        self.name = 'AimsCube'
        #self.origin = origin
        #self.edges = edges
        #self.points = points
        self.plots = plots 

    def ncubes(self):
        """returns the number of cube files to output """
        if self.plots:
            number = len(self.plots)
        else:
            number = 0
        return number

    def set(self, **kwargs):
        """ set any of the parameters ... """
        # NOT IMPLEMENTED AT THE MOMENT!

    #def move_to_base_name(self, basename):
        #""" when output tracking is on or the base namem is not standard,
        #this routine will rename add the base to the cube file output for
        #easier tracking """
        #for plot in self.plots:
            #found = False
            #cube = plot.split()
            #if (cube[0] == 'total_density' or
                #cube[0] == 'spin_density' or
                #cube[0] == 'delta_density'):
                #found = True
                #old_name = cube[0] + '.cube'
                #new_name = basename + '.' + old_name
            #if cube[0] == 'eigenstate' or cube[0] == 'eigenstate_density':
                #found = True
                #state = int(cube[1])
                #s_state = cube[1]
                #for i in [10, 100, 1000, 10000]:
                    #if state < i:
                        #s_state = '0' + s_state
                #old_name = cube[0] + '_' + s_state + '_spin_1.cube'
                #new_name = basename + '.' + old_name
            #if found:
                #os.system('mv ' + old_name + ' ' + new_name)
#
    def add_plot(self, name):
        """ in case you forgot one ... """
        self.plots += [name]

    def write(self, file):
        """ write the necessary output to the already opened control.in """
        file.write('output cube ' + self.plots[0] + '\n')
        file.write('output cube ' + self.plots[0] + '\n')
        file.write('cube spinstate 2\n')
        #file.write('   cube origin ')
        #for ival in self.origin:
            #file.write(str(ival) + ' ')
        #file.write('\n')
        #for i in range(3):
            #file.write('   cube edge ' + str(self.points[i]) + ' ')
            #for ival in self.edges[i]:
                #file.write(str(ival) + ' ')
            #file.write('\n')
        if self.ncubes() > 1:
            for i in range(self.ncubes() - 1):
                file.write('output cube ' + self.plots[i + 1] + '\n')
                file.write('output cube ' + self.plots[i + 1] + '\n')
                file.write('cube spinstate 2\n')

#class fo_cpmd(fodft):
#    """ Wrapper for old cpmd io/calculator. Need to integrate into this later on..."""
#   
#    def __init__(self, dimer, fformat="xyz", charges=1):
#        fodft.__init__(self, dimer, fformat="xyz")
#
#        self.dimer.calc = CPMD(self.dimer)
#        self.dimer.calc.center(vacuum=4)
#        self.dimer.calc.set(charge=charges)
#
#    def create_fragments(self, fo_files=True):
#        self.dimer.calc.create_fragments(self.dimer, fo_files=True)


