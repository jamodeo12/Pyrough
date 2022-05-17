# ---------------------------------------------------------------------------
# 
# Title: Main code pyrough
#
# Authors: Jonathan Amodeo Javier Gonzalez Jennifer Izaguirre Christophe Le Bourlot
#
# Date: January 30, 2020
#
# Pyrough is a code used to provide the user with either FEM or MD files of various objects that have had random surface roughness applied to it.
# The objects included in this code are nanowire , slab, sphere. In the main code the user must provide what type of 
# object thy desired. Based on this the code will read the parameters from the json file provided and return an stl if FEM is desired or
# lmp file if MD is desired.  
# 
# ---------------------------------------------------------------------------
from Sources import Param_class
from Sources import Sample_class
import sys, subprocess

print('#########################################################################')
print('#                    Running Random Surface script                      #')
print('#     Jonathan Amodeo & Javier Gonzalez & Jennifer Izaguirre 2020       #')
print('#########################################################################')

# _____________________Removing previous data____________________

subprocess.call(['rm', 'material_supercell.lmp'])
subprocess.call(['rm', 'sample_with_atoms.lmp'])

# _____________________Main Code____________________

Param_file = sys.argv[1]

param = Param_class.Parameter(Param_file)

# -----------------------second task callingthe nanowire into the sample to get an stl

# SAMPLE IS THE PATH GUIDE TO CALLING THE CLASS TO CREATE THE STL FILE

sample = Sample_class.Sample(param.type_S)  # THIS CLASS CALLING REAUIRES A SAMPLE WITH A DESIRED OBJECT NAME AS AN INPUT

vertices, FEM_stl = sample.make_stl(param.type_S,
                          param.B,
                          param.C1,
                          param.N,
                          param.M,
                          param.radius,
                          param.length,
                          param.height,
                          param.width,
                          param.ns,
                          param.raw_stl,
                          param.nfaces,
                          param.surfaces,
                          param.energies,
                          param.n_at,
                          param.lattice_structure,
                          param.lattice_parameter,
                          param.material)
# Calling the sample class function of MAKE it which returns an stl file of the object desired.
if param.output(Param_file) == 'ATOM_lmp':
    sample.make_atom(FEM_stl,
                     param.lattice_structure,
                     param.lattice_parameter,
                     param.material,
                     param.orien_x,
                     param.orien_y,
                     param.orien_z,
                     vertices)
    # call make it md to create atomsk file
    print('JOB DONE!' + '  File name: sample_with_atoms.lmp')
else:
    print('JOB DONE!' + '  File name: ' + FEM_stl)
