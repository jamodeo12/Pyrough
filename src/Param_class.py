# ---------------------------------------------------------------------------
# Title: PARAM CLASS
# Authors: Jonathan Amodeo, Hugo Iteney, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot
# Date: June 01, 2022
#
# The PARAM class is responsible for reading/storing user parameters from a JSON file.
# The class also reads the desired file format requested by the user, either for atomistic or FEM exports.
# IMPORTANT: If any values need to be modified, update the JSON file directly.
# Do not modify the variable storage within the code.
# ---------------------------------------------------------------------------

import json
from src import Func_pyrough as fp

class Parameter:
    """
    This class stores the parameters defined by the user using the JSON input file.

    """

    def __init__(self, json_file):
        with open(json_file) as json_file:
            read_param = json.load(json_file)

        if "Precinmatrix" in read_param:
            self.spec = read_param["Precinmatrix"]["Spec"]
            self.length_x2 = read_param["MATRIX_Param"]["Length_x"]
            self.length_y2 = read_param["MATRIX_Param"]["Length_y"]
            self.length_z2 = read_param["MATRIX_Param"]["Length_z"]
            self.lattice_structure2 = read_param["MATRIX_Param"]["Lattice_structure"]
            #self.lattice_structure2 = fp.convert_in_list_of_string(self.lattice_structure2)
            self.lattice_parameter2 = read_param["MATRIX_Param"]["Lattice_parameter"]
            self.lattice_parameter2 = fp.convert_in_list_of_string(self.lattice_parameter2)
            self.material2 = read_param["MATRIX_Param"]["Material"]
            self.material2= fp.convert_in_list_of_string(self.material2)
            self.orien_x2 = read_param["MATRIX_Param"]["Orien_x"]
            self.orien_y2 = read_param["MATRIX_Param"]["Orien_y"]
            self.orien_z2 = read_param["MATRIX_Param"]["Orien_z"]
            try:
                self.precpos = read_param["MATRIX_Param"]["Precpos"]
            except KeyError:
                self.precpos = 'center'

        if "Grain" in read_param:
            self.type_S = "grain"
            self.N = read_param["Grain"]["N"]
            self.M = read_param["Grain"]["M"]
            self.eta = read_param["Grain"]["eta"]
            try:
                self.C1 = read_param["Grain"]["C1"]
            except KeyError:
                self.C1 = ""
            try:
                self.RMS = read_param["Grain"]["RMS"]
            except KeyError:
                self.RMS = ""
            self.height = read_param["Grain"]["Height"]
            self.length = read_param["Grain"]["Length"]
            self.width = read_param["Grain"]["Width"]
            self.ns = read_param["Grain"]["Mesh_size"]
            try:
                self.alpha = read_param["Grain"]["Refine_factor"]
            except KeyError:
                self.alpha = 1
            self.raw_stl = read_param["Grain"]["Raw_stl"]
            self.lattice_structure1 = read_param["ATOM1_Param"]["Lattice_structure"]
            self.lattice_parameter1 = read_param["ATOM1_Param"]["Lattice_parameter"]
            self.lattice_parameter1 = fp.convert_in_list_of_string(self.lattice_parameter1)
            self.material1 = read_param["ATOM1_Param"]["Material"]
            self.material1 = fp.convert_in_list_of_string(self.material1)
            self.orien_x1 = read_param["ATOM1_Param"]["Orien_x"]
            self.orien_y1 = read_param["ATOM1_Param"]["Orien_y"]
            self.orien_z1 = read_param["ATOM1_Param"]["Orien_z"]
            self.lattice_structure2 = read_param["ATOM2_Param"]["Lattice_structure"]
            self.lattice_parameter2 = read_param["ATOM2_Param"]["Lattice_parameter"]
            self.lattice_parameter2 = fp.convert_in_list_of_string(self.lattice_parameter2)
            self.material2 = read_param["ATOM2_Param"]["Material"]
            self.material2 = fp.convert_in_list_of_string(self.material2)
            self.orien_x2 = read_param["ATOM2_Param"]["Orien_x"]
            self.orien_y2 = read_param["ATOM2_Param"]["Orien_y"]
            self.orien_z2 = read_param["ATOM2_Param"]["Orien_z"]
            self.ext_fem = read_param["Output"]["FEM"]
            self.ext_ato = read_param["Output"]["ATOM"]

        elif "Multi_layered" in read_param:
            self.type_S = "multi_layered"
            self.N = read_param["Multi_layered"]["N"]
            self.M = read_param["Multi_layered"]["M"]
            self.eta = read_param["Multi_layered"]["eta"]
            try:
                self.C1 = read_param["Multi_layered"]["C1"]
            except KeyError:
                self.C1 = ""
            try:
                self.RMS = read_param["Multi_layered"]["RMS"]
            except KeyError:
                self.RMS = ""
            self.height = read_param["Multi_layered"]["Height"]
            self.length = read_param["Multi_layered"]["Length"]
            self.width = read_param["Multi_layered"]["Width"]
            self.ns = read_param["Multi_layered"]["Mesh_size"]
            try:
                self.alpha = read_param["Multi_layered"]["Refine_factor"]
            except KeyError:
                self.alpha = 1
            self.raw_stl = read_param["Multi_layered"]["Raw_stl"]
            self.lattice_structure = read_param["ATOMS_Param"]["Lattice_structure"]
            self.lattice_parameter = read_param["ATOMS_Param"]["Lattice_parameter"]
            for lp in range(len(self.lattice_parameter)):
                self.lattice_parameter[lp] = fp.convert_in_list_of_string(self.lattice_parameter[lp])
            self.material = read_param["ATOMS_Param"]["Material"]
            for mat in range(len(self.material)):
                self.material[mat] = fp.convert_in_list_of_string(self.material[mat])
            self.orien_x = read_param["ATOMS_Param"]["Orien_x"]
            self.orien_y = read_param["ATOMS_Param"]["Orien_y"]
            self.orien_z = read_param["ATOMS_Param"]["Orien_z"]
            try:
                self.pattern_layer = read_param["Multi_layered"]["pattern_layer"]
            except KeyError:
                self.pattern_layer = []
                for i in range(len(self.material)):
                    self.pattern_layer.append(i)
            try:
                self.height_layer = read_param["Multi_layered"]["height_layer"]
            except KeyError:
                self.height_layer = []
                j=0
                while j < (len(self.material)):
                    self.height_layer.append(self.height/len(self.pattern_layer))
                    j=j+1

            self.ext_fem = read_param["Output"]["FEM"]
            self.ext_ato = read_param["Output"]["ATOM"]

        else:
            if "cWire" in read_param:
                self.type_S = "cwire"
                self.eta = read_param["cWire"]["eta"]
                try:
                    self.C1 = read_param["cWire"]["C1"]
                except KeyError:
                    self.C1 = ""
                try:
                    self.RMS = read_param["cWire"]["RMS"]
                except KeyError:
                    self.RMS = ""
                self.N = read_param["cWire"]["N"]
                self.M = read_param["cWire"]["M"]
                self.length = read_param["cWire"]["Length"]
                self.radius = read_param["cWire"]["Radius"]
                self.ns = read_param["cWire"]["Mesh_size"]
                try:
                    self.alpha = read_param["cWire"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param["cWire"]["Raw_stl"]
                try:
                    self.angles = read_param["cWire"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "Box" in read_param:
                self.type_S = "box"
                self.N = read_param["Box"]["N"]
                self.M = read_param["Box"]["M"]
                self.eta = read_param["Box"]["eta"]
                try:
                    self.C1 = read_param["Box"]["C1"]
                except KeyError:
                    self.C1 = ""
                try:
                    self.RMS = read_param["Box"]["RMS"]
                except KeyError:
                    self.RMS = ""
                self.height = read_param["Box"]["Height"]
                self.length = read_param["Box"]["Length"]
                self.width = read_param["Box"]["Width"]
                self.ns = read_param["Box"]["Mesh_size"]
                try:
                    self.alpha = read_param["Box"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.radius = 0
                self.raw_stl = read_param["Box"]["Raw_stl"]
                try:
                    self.angles = read_param["Box"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "Sphere" in read_param:
                self.type_S = "sphere"
                self.eta = read_param["Sphere"]["eta"]
                self.N = read_param["Sphere"]["N"]
                self.M = read_param["Sphere"]["M"]
                self.C1 = read_param["Sphere"]["C1"]
                self.RMS = ""
                self.radius = read_param["Sphere"]["Radius"]
                self.ns = read_param["Sphere"]["Mesh_size"]
                try:
                    self.alpha = read_param["Sphere"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.length = 0
                self.width = 0
                self.raw_stl = read_param["Sphere"]["Raw_stl"]
                try:
                    self.angles = read_param["Sphere"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "fWire" in read_param:
                self.type_S = "fwire"
                self.eta = read_param["fWire"]["eta"]
                try:
                    self.C1 = read_param["fWire"]["C1"]
                except KeyError:
                    self.C1 = ""
                try:
                    self.RMS = read_param["fWire"]["RMS"]
                except KeyError:
                    self.RMS = ""
                self.N = read_param["fWire"]["N"]
                self.M = read_param["fWire"]["M"]
                self.length = read_param["fWire"]["Length"]
                self.nfaces = read_param["fWire"]["N_Faces"]
                self.radius = read_param["fWire"]["Radius"]
                self.ns = read_param["fWire"]["Mesh_size"]
                try:
                    self.alpha = read_param["fWire"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param["fWire"]["Raw_stl"]
                try:
                    self.angles = read_param["fWire"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "Wulff" in read_param:
                self.type_S = "wulff"
                self.eta = read_param["Wulff"]["eta"]
                try:
                    self.C1 = read_param["Wulff"]["C1"]
                except KeyError:
                    self.C1 = ""
                try:
                    self.RMS = read_param["Wulff"]["RMS"]
                except KeyError:
                    self.RMS = ""
                self.n_at = read_param["Wulff"]["n_atoms"]
                self.surfaces = read_param["Wulff"]["Surfaces"]
                self.energies = read_param["Wulff"]["Energies"]
                self.N = read_param["Wulff"]["N"]
                self.M = read_param["Wulff"]["M"]
                self.length = 0
                self.nfaces = 0
                self.radius = 0
                self.ns = read_param["Wulff"]["Mesh_size"]
                try:
                    self.alpha = read_param["Wulff"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param["Wulff"]["Raw_stl"]
                try:
                    self.angles = read_param["Wulff"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "Cube" in read_param:
                self.type_S = "cube"
                self.eta = read_param["Cube"]["eta"]
                try:
                    self.C1 = read_param["Cube"]["C1"]
                except KeyError:
                    self.C1 = ""
                try:
                    self.RMS = read_param["Cube"]["RMS"]
                except KeyError:
                    self.RMS = ""
                self.n_at = 0
                self.surfaces = 0
                self.energies = 0
                self.N = read_param["Cube"]["N"]
                self.M = read_param["Cube"]["M"]
                self.length = read_param["Cube"]["Length"]
                self.nfaces = 0
                self.radius = 0
                self.ns = read_param["Cube"]["Mesh_size"]
                try:
                    self.alpha = read_param["Cube"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = read_param["Cube"]["Length"]
                self.width = read_param["Cube"]["Length"]
                self.raw_stl = read_param["Cube"]["Raw_stl"]
                try:
                    self.angles = read_param["Cube"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                #self.angles = read_param["Cube"]["Angles"]
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

#------------------------------------------------------------------------------------------------------------------------------
### Pillar part: under development
            elif "fPillar" in read_param:
                self.type_S = "fpillar"
                self.eta = read_param["fPillar"]["eta"]
                try:
                    self.C1 = read_param["fPillar"]["C1"]
                except KeyError:
                    self.C1 =["",""]
                try:
                    self.RMS = read_param["fPillar"]["RMS"]
                except KeyError:
                    self.RMS =["",""]
                self.N = read_param["fPillar"]["N"]
                self.M = read_param["fPillar"]["M"]
                self.length = read_param["fPillar"]["Length"]
                self.nfaces = read_param["fPillar"]["N_Faces"]
                self.radius = read_param["fPillar"]["Radius"]
                self.ns = read_param["fPillar"]["Mesh_size"]
                try:
                    self.alpha = read_param["fPillar"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param["fPillar"]["Raw_stl"]
                try:
                    self.angles = read_param["fPillar"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

            elif "cPillar" in read_param:
                self.type_S = "cpillar"
                self.eta = read_param["cPillar"]["eta"]
                try:
                    self.C1 = read_param["cPillar"]["C1"]
                except KeyError:
                    self.C1 =["",""]
                try:
                    self.RMS = read_param["cPillar"]["RMS"]
                except KeyError:
                    self.RMS =["",""]
                self.N = read_param["cPillar"]["N"]
                self.M = read_param["cPillar"]["M"]
                self.length = read_param["cPillar"]["Length"]
                self.nfaces = read_param["cPillar"]["N_Faces"]
                self.radius = read_param["cPillar"]["Radius"]
                self.ns = read_param["cPillar"]["Mesh_size"]
                try:
                    self.alpha = read_param["cPillar"]["Refine_factor"]
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param["cPillar"]["Raw_stl"]
                try:
                    self.angles = read_param["cPillar"]["Angles"]
                except KeyError:
                    self.angles = [0,0,0]
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param["Output"]["FEM"]
                self.ext_ato = read_param["Output"]["ATOM"]

#-------------------------------------------------------------------------------------------------------------------------------------------------
            self.lattice_structure = read_param["ATOM_Param"]["Lattice_structure"]
            self.lattice_parameter = read_param["ATOM_Param"]["Lattice_parameter"]
            self.lattice_parameter = fp.convert_in_list_of_string(self.lattice_parameter)
            self.material = read_param["ATOM_Param"]["Material"]
            self.material = fp.convert_in_list_of_string(self.material)
            self.orien_x = read_param["ATOM_Param"]["Orien_x"]
            self.orien_y = read_param["ATOM_Param"]["Orien_y"]
            self.orien_z = read_param["ATOM_Param"]["Orien_z"]
            try:
                self.angles2 = read_param["ATOM_Param"]["Angles"]
            except KeyError:
                self.angles2 = [0, 0, 0]

    def output(self, json_file):
        """
        This function obtains the desired file type of the user. Having the option to obtain an stl
        file or lmp file based on the decision if MD or FEM is desired.

        :param json_file: A json file that that has the parameters stored into it
        :type json_file: file
        """
        with open(json_file) as f:
            read_param = json.load(f)
        if len(read_param["Output"]["ATOM"]) > 0:
            return "ATOM"
