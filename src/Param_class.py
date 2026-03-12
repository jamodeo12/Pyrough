# ---------------------------------------------------------------------------
# Title: PARAM CLASS
# Authors: Jonathan Amodeo, Hugo Iteney, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot, Thomas Cadart
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

        if "Output" in read_param:
            self.ext_fem = read_param["Output"].get("FEM", "")
            self.ext_ato = read_param["Output"].get("ATOM", "")

        if "ATOM_Param" in read_param:
            atom_param = read_param["ATOM_Param"]

            def to_list(val):
                """Normalize to list: scalar → [scalar], single item → [item]"""
                return val if isinstance(val, list) else [val]

            def to_list_of_lists(val):
                """Normalize orientations/parameters to always be a list of lists.
                [1,0,0]        → [[1,0,0]]        (single orientation)
                [[1,0,0],[...]]→ [[1,0,0],[...]]  (multiple orientations, unchanged)
                3.615          → [[3.615]]        (single scalar parameter)
                [3.615, 4.07]  → [[3.615],[4.07]] (multiple scalar parameters)
                """
                if not isinstance(val, list):
                    return [[val]]
                if len(val) == 0:
                    return [val]
                # If first element is a list → already list of lists
                if isinstance(val[0], list):
                    return val
                # If first element is a number and length matches an orientation (3 or 4)
                if isinstance(val[0], (int, float)) and len(val) in (3, 4):
                    return [val]  # single orientation vector
                # Otherwise it's a list of scalar parameters → wrap each
                return [[v] if not isinstance(v, list) else v for v in val]

            self.lattice_structure = to_list(atom_param["Lattice_structure"])
            self.lattice_parameter = [fp.convert_in_list_of_string(lp)
                                      for lp in to_list_of_lists(atom_param["Lattice_parameter"])]
            self.material = [fp.convert_in_list_of_string(m)
                             for m in to_list_of_lists(atom_param["Material"])]
            self.orien_x = to_list_of_lists(atom_param["Orien_x"])
            self.orien_y = to_list_of_lists(atom_param["Orien_y"])
            self.orien_z = to_list_of_lists(atom_param["Orien_z"])
            self.angles2 = atom_param.get("Angles", [0, 0, 0])

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
                self.precpos = read_param["MATRIX_Param"].get("Precpos", 'center')

        if self.type_S == "lattice" and "Output" in read_param:
            try:
                self.geometry_lattice = read_param["geometry_lattice"]
            except KeyError:
                raise ValueError("No geometry_lattice key found in the JSON file for lattice sample type.")
        if self.type_S == "lattice":
            self.trim_boundary_beams = read_param.get("trim_boundary_beams", False)
            if self.trim_boundary_beams:
                self.cube_trimmer_param = read_param.get("cube_trimmer_parameters", None)
                if self.cube_trimmer_param is None:
                    print("No cube_trimmer_parameters key found in the JSON file for lattice sample type. The Wire "
                          "parameters will be used as default values for the cube trimmer.")

        # Check that all required parameters are present for the given sample type
        if "Output" in read_param:
            self._check_required_params()

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

    def _check_required_params(self):
        """
        Checks that all required parameters are present for the given sample type.

        :raises ValueError: If any required parameter is missing or has a default/empty value
        """

        required_sample_params = {
            "box": ["N", "M", "eta", "height", "width", "length", "ns", "alpha"],
            "cube": ["N", "M", "eta", "height", "width", "length", "ns", "alpha"],
            "cwire": ["N", "M", "eta", "length", "radius", "ns", "alpha"],
            "fwire": ["N", "M", "eta", "length", "radius", "nfaces", "ns", "alpha"],
            "fpillar": ["N", "M", "eta", "length", "radius", "nfaces", "ns", "alpha"],
            "cpillar": ["N", "M", "eta", "length", "radius", "ns", "alpha"],
            "pillar": ["N", "M", "eta", "length", "radius", "ns", "alpha"],
            "sphere": ["N", "M", "eta", "radius", "ns", "alpha"],
            "wulff": ["N", "M", "eta", "n_at", "surfaces", "energies", "ns", "alpha"],
            "grain": ["N", "M", "eta", "height", "width", "length", "ns", "alpha"],
            "multi_layered": ["N", "M", "eta", "height", "width", "length", "ns", "alpha",
                              "pattern_layer", "height_layer"],
            "lattice": ["N", "M", "eta", "ns", "alpha", "beam_type"],
        }

        required_atom_params = ["lattice_structure", "lattice_parameter", "material",
                                "orien_x", "orien_y", "orien_z"]

        required_lattice_params = ["geometry_lattice"]

        errors = []

        # Check sample parameters
        missing_sample = [p for p in required_sample_params.get(self.type_S, [])
                          if not getattr(self, p, None)]
        if missing_sample:
            errors.append(f"  - '{self.key}' block: {', '.join(missing_sample)}")

        # Check atom parameters
        missing_atom = [p for p in required_atom_params if not getattr(self, p, None)]
        if missing_atom:
            errors.append(f"  - 'ATOM_Param' block: {', '.join(missing_atom)}")

        # Check lattice-specific parameters
        if self.type_S == "lattice":
            missing_lattice = [p for p in required_lattice_params if not getattr(self, p, None)]
            if missing_lattice:
                errors.append(f"  - 'geometry_lattice' block: {', '.join(missing_lattice)}")

        if errors:
            raise ValueError("Missing parameters in JSON file:\n" + "\n".join(errors))

