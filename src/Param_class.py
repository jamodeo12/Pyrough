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

        type_allowed = ["grain", "multi_layered", "cwire", "fwire", "box", "sphere", "poly", "wulff", "cube", "fpillar",
                        "pillar", "poly_pillar", "lattice"]

        if not any(key.lower() in type_allowed for key in read_param):
            raise ValueError("No valid sample type found in the JSON file.")
        else:
            for key in read_param:
                if key.lower() in type_allowed:
                    self.type_S = key.lower()
                    self.key = key

        self.N = read_param[self.key].get("N", 0)
        self.M = read_param[self.key].get("M", 0)
        self.eta = read_param[self.key].get("eta", 0.0)
        if type(self.eta) is int or type(self.eta) is float:
            if self.type_S == "sphere":
                self.B = 4 * self.eta
            else:
                self.B = 2 * (1 + self.eta)
        self.C1 = read_param[self.key].get("C1", "")
        self.RMS = read_param[self.key].get("RMS", "")
        self.radius = read_param[self.key].get("Radius", 0.0)
        self.height = read_param[self.key].get("Height", 0.0)
        self.length = read_param[self.key].get("Length", 0.0)
        self.width = read_param[self.key].get("Width", 0.0)
        self.nfaces = read_param[self.key].get("N_Faces", 0)
        self.surfaces = read_param[self.key].get("Surfaces", [])
        self.energies = read_param[self.key].get("Energies", [])
        self.n_at = read_param[self.key].get("n_atoms", 0)
        self.ns = read_param[self.key].get("Mesh_size", 0.0)
        self.angles = read_param[self.key].get("Angles", [0, 0, 0])
        self.alpha = read_param[self.key].get("Refine_factor", 1)
        self.beam_type = read_param[self.key].get("beam_type", "")
        self.raw_stl = read_param[self.key].get("Raw_stl", "")
        self.ext_fem = read_param["Output"]["FEM"]
        self.ext_ato = read_param["Output"]["ATOM"]

        if "ATOM_Param" in read_param:
            self.lattice_structure = read_param["ATOM_Param"]["Lattice_structure"]
            self.lattice_parameter = read_param["ATOM_Param"]["Lattice_parameter"]
            self.lattice_parameter = fp.convert_in_list_of_string(self.lattice_parameter)
            self.material = read_param["ATOM_Param"]["Material"]
            self.material = fp.convert_in_list_of_string(self.material)
            self.orien_x = read_param["ATOM_Param"]["Orien_x"]
            self.orien_y = read_param["ATOM_Param"]["Orien_y"]
            self.orien_z = read_param["ATOM_Param"]["Orien_z"]
        elif "ATOM1_Param" in read_param and "ATOM2_Param" in read_param:
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
        elif "ATOMS_Param" in read_param:
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
            self.angles2 = read_param["ATOMS_Param"].get("Angles", [0, 0, 0])

        if "Multi_layered" in read_param:
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
                j = 0
                while j < (len(self.material)):
                    self.height_layer.append(self.height / len(self.pattern_layer))
                    j += 1

            if self.type_S == "poly_pillar" or self.type_S == "pillar":
                self.C1 = ["", ""]
                self.RMS = ["", ""]

            if "Precinmatrix" in read_param:
                self.spec = read_param["Precinmatrix"]["Spec"]
                self.length_x2 = read_param["MATRIX_Param"]["Length_x"]
                self.length_y2 = read_param["MATRIX_Param"]["Length_y"]
                self.length_z2 = read_param["MATRIX_Param"]["Length_z"]
                self.lattice_structure2 = read_param["MATRIX_Param"]["Lattice_structure"]
                # self.lattice_structure2 = fp.convert_in_list_of_string(self.lattice_structure2)
                self.lattice_parameter2 = read_param["MATRIX_Param"]["Lattice_parameter"]
                self.lattice_parameter2 = fp.convert_in_list_of_string(self.lattice_parameter2)
                self.material2 = read_param["MATRIX_Param"]["Material"]
                self.material2 = fp.convert_in_list_of_string(self.material2)
                self.orien_x2 = read_param["MATRIX_Param"]["Orien_x"]
                self.orien_y2 = read_param["MATRIX_Param"]["Orien_y"]
                self.orien_z2 = read_param["MATRIX_Param"]["Orien_z"]

        if self.type_S == "lattice":
            try:
                self.geometry_lattice = read_param["geometry_lattice"]
            except KeyError:
                raise ValueError("No geometry_lattice key found in the JSON file for lattice sample type.")
            self.trim_boundary_beams = read_param["geometry_lattice"].get("trim_boundary_beams", False)
            if self.trim_boundary_beams:
                self.cube_trimmer_param = read_param["geometry_lattice"].get("cube_trimmer_parameters", None)
                if self.cube_trimmer_param is None:
                    print("No cube_trimmer_parameters key found in the JSON file for lattice sample type. The Wire "
                          "parameters will be used as default values for the cube trimmer.")


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
