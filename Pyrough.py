# ---------------------------------------------------------------------------
# 
# Title: Main code pyrough
#
# Authors: Jonathan Amodeo, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot & Hugo Iteney
#
# Date: June 01, 2022
#
# Pyrough is a code used to provide the user with either FEM or MD files of various objects that have had random surface roughness applied to it.
# The objects included in this code are nanowire , slab, sphere. In the main code the user must provide what type of 
# object thy desired. Based on this the code will read the parameters from the json file provided and return an stl if FEM is desired or
# lmp file if MD is desired.  
# 
# ---------------------------------------------------------------------------
from src import Param_class
from src import Sample_class
import os
import sys, subprocess

print('#######################################################################################')
print('#                                       Pyrough                                       #')
print('#                          Jonathan Amodeo & Hugo Iteney 2023                         #')
print('#######################################################################################')

# _____________________Removing previous data____________________

# subprocess.call(['rm', 'material_supercell.lmp'])
# subprocess.call(['rm', 'sample_with_atoms.lmp'])

# _____________________Main Code____________________
if sys.argv[1] == '-surface' :
    current_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess.call(['python', current_dir+'/src/Surface_Analysis.py', sys.argv[2], sys.argv[3], sys.argv[4]])
else :
    Param_file = sys.argv[1]
    out_pre = os.path.basename(Param_file)[:-5]

    param = Param_class.Parameter(Param_file)

    # -----------------------second task calling the nanowire into the sample to get an stl

    # SAMPLE IS THE PATH GUIDE TO CALLING THE CLASS TO CREATE THE STL FILE

    if param.type_S == 'grain' :
        vertices, FEM_stl = Sample_class.make_box(param.type_S,
                                                  2*(1+param.H),
                                                  param.C1,
                                                  param.RMS,
                                                  param.N,
                                                  param.M,
                                                  param.length,
                                                  param.height,
                                                  param.width,
                                                  param.ns,
                                                  param.raw_stl,
                                                  out_pre)
        Sample_class.make_atom_grain(FEM_stl,
                             param.lattice_structure1,
                             param.lattice_parameter1,
                             param.material1,
                             param.orien_x1,
                             param.orien_y1,
                             param.orien_z1,
                             param.lattice_structure2,
                             param.lattice_parameter2,
                             param.material2,
                             param.orien_x2,
                             param.orien_y2,
                             param.orien_z2,
                             vertices,
                             out_pre)
            # call make it md to create atomsk file
        print('JOB DONE!' + '  File name: ' + out_pre + '.lmp')

    else :
        sample = Sample_class.Sample(param.type_S)  # THIS CLASS CALLING REAUIRES A SAMPLE WITH A DESIRED OBJECT NAME AS AN INPUT

        vertices, FEM_stl = sample.make_stl(param.type_S,
                                            param.H,
                                            param.C1,
                                            param.RMS,
                                            param.N,
                                            param.M,
                                            param.radius,
                                            param.length,
                                            param.height,
                                            param.width,
                                            param.ns,
                                            param.alpha,
                                            param.raw_stl,
                                            param.nfaces,
                                            param.surfaces,
                                            param.energies,
                                            param.n_at,
                                            param.lattice_structure,
                                            param.lattice_parameter,
                                            param.material,
                                            param.orien_x,
                                            param.orien_y,
                                            param.orien_z,
                                            out_pre,
                                            param.ext_fem)
        print('====== > FEM JOB DONE !')
        for ext in param.ext_fem:
            print('====== > File name: ' + out_pre + '.' + str(ext))
        if 'stl' not in param.ext_fem :
            print('====== > File name: ' + out_pre + '.' + str(ext))
        # Calling the sample class function of MAKE it which returns an stl file of the object desired.
        if param.output(Param_file) == 'ATOM':
            sample.make_atom(param.type_S,
                             FEM_stl,
                             param.lattice_structure,
                             param.lattice_parameter,
                             param.material,
                             param.orien_x,
                             param.orien_y,
                             param.orien_z,
                             vertices,
                             out_pre,
                             param.ext_ato)
            # call make it md to create atomsk file
            # print('JOB DONE!' + '  File name: ' + out_pre + '.' + str(param.ext_ato[0]))
            print('====== > ATOM JOB DONE !')
            for ext in param.ext_ato:
                print('====== > File name: ' + out_pre + '.' + str(ext))
        # else:
            # print('JOB DONE!' + '  File name: ' + FEM_stl)
