# ---------------------------------------------------------------------------
# Title: Main code Pyrough
# Authors: Jonathan Amodeo, Hugo Iteney, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot
# Date: June 01, 2022
#
# Pyrough generates 3D virtual samples with rough surfaces for atomistic and finite-element simulations.
# Several objects such as spheres, cubes, Wulff-shaped systems or wires can be generated.
# Additional options allowing to generate rough GBs or post-treat experimental images are also available.
# Various output formats are handled including .stl, .msh, .xyz, .lmp, etc.
# Example input .json files are provided in the examples/ folder, have a look !
# ---------------------------------------------------------------------------
import os
import subprocess
import sys
from src import Param_class, Sample_class, Func_pyrough as fp

print("##################################################################")
print("#    /\/\/\                                                      #")
print("#    | O _ |        Pyrough                                      #")
print("#  <|-|--/          Version 1.2                                  #")
print("#      \---|-|>     (C) 2023 Jonathan Amodeo and co.             #")
print("#      | _ O |      https://github.com/jamodeo12/Pyrough         #")
print("#      \/\/\/                                                    #")
print("##################################################################")

# _____________________Main Code____________________
if sys.argv[1] == "-surface":
    # -surface : Image characterization and digital twin generation
    print("====== > Pyrough.py : surface analysis treatment running...")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess.call(
        [
            "python",
            current_dir + "/src/Surface_Analysis.py",
            sys.argv[2],
            sys.argv[3],
            sys.argv[4],
            sys.argv[5],
        ]
    )
elif sys.argv[1] == "-test_pyrough_execution":
    # -test_pyrough_execution : test if all examples/*.json files can be processed
    print("====== > Pyrough.py : test_pyrough_execution running...")
    print("====== > Pyrough.py : check/execute all .json files in examples folder")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    fp.test_pyrough_execution(current_dir)
else:
    # If no specific option, we generate a rough sample
    # Parameter file definition
    print("====== > Pyrough.py : parameter file definition running...")
    Param_file = sys.argv[1]
    out_pre = os.path.basename(Param_file)[:-5]

    param = Param_class.Parameter(Param_file)

    if param.type_S == "grain":
        # If grain, we first create the mesh of a rough box
        # the STL file will be used to make grain 1, while its negative will be used for grain 2
        print("====== > Pyrough.py : grain option treatment running...")
        print("====== > Pyrough.py : Sample_class.make_box running...")
        vertices, FEM_stl = Sample_class.make_box(
            param.type_S,
            2 * (1 + param.eta),
            param.C1,
            param.RMS,
            param.N,
            param.M,
            param.length,
            param.height,
            param.width,
            param.ns,
            param.alpha,
            param.raw_stl,
            out_pre,
            param.ext_fem,
        )

        print("====== > Pyrough.py : Sample_class.make_atom_grain running...")
        Sample_class.make_atom_grain(
            FEM_stl,
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
            out_pre,
            param.ext_ato,
        )
        # call make it md to create atomsk file
        print("JOB DONE!" + "  File name: " + out_pre + ".lmp")

    else:
        print("====== > Pyrough.py : {} option treatment running...".format(param.type_S))
        sample = Sample_class.Sample(param.type_S)

        print("====== > Pyrough.py : sample.make_stl running...")
        vertices, FEM_stl = sample.make_stl(
            param.type_S,
            param.eta,
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
            param.orien_z,
            out_pre,
            param.ext_fem,
        )

        print("====== > FEM JOB DONE !")
        for ext in param.ext_fem:
            print("====== > File name: " + out_pre + "." + str(ext))

        if "stl" not in param.ext_fem:
            print("====== > File name: " + out_pre + ".stl")

        if param.output(Param_file) == "ATOM" and not hasattr(param, "spec"):
            print("====== > Pyrough.py : sample.make_atom running...")
            sample.make_atom(
                FEM_stl,
                param.lattice_structure,
                param.lattice_parameter,
                param.material,
                param.orien_x,
                param.orien_y,
                param.orien_z,
                vertices,
                out_pre,
                param.ext_ato,
            )

            # call make it md to create atomsk file
            print("====== > ATOM OBJECT JOB DONE !")
            for ext in param.ext_ato:
                print("====== > File name: " + out_pre + "." + str(ext))

        if  hasattr(param, "spec"):
             if param.spec == "precinmatrix":
                print("====== > Pyrough.py : precinmatric option treatment running...")

                print("====== > sample.make_precipitate running")
                # Generate the precipitate centered in the matrix-sized supercell : precipitate.lmp
                sample.make_precipitate(
                    FEM_stl,
                    param.lattice_structure,
                    param.lattice_parameter,
                    param.material,
                    param.orien_x,
                    param.orien_y,
                    param.orien_z,
                    param.length_x2,
                    param.length_y2,
                    param.length_z2,
                )

                print("====== > sample.make_atom_matrix running")
                sample.make_atom_matrix(
                FEM_stl,
                param.length_x2,
                param.length_y2,
                param.length_z2,
                param.lattice_structure2,
                param.lattice_parameter2,
                param.material2,
                param.orien_x2,
                param.orien_y2,
                param.orien_z2,
                )

                print("====== > sample.put_prec_in_matrix running...")
                sample.put_prec_in_matrix(
                    out_pre,
                    param.ext_ato
                 )

                print("====== > ATOM PRECIPITATE IN MATRIX JOB DONE !")
