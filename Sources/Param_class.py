# ---------------------------------------------------------------------------
# 
# Title: PARAM CLASS
#
# Authors: Jonathan Amodeo, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot & Hugo Iteney
#
# Date: June 01, 2022
#
# This class file is used to store the parameters desired by the user. 
# It requires a JSON file to be able to function, as it reads directly
# from the file and stores the values in a variable.
# The code also reads the desired file format the user is requesting, either FEM or MD.
# IMPORTANT! If values need to be changed refer to the JSON file. Do not change 
# the storage of the variable.
# ---------------------------------------------------------------------------
# This command opens the json file and reads what the user has inputted
# The user needs to specify the type of sample they desire in the _int_ for this loop to wor
# Variable Declaration based on parameters in the json file


import json


class Parameter(object):
    """
    This class is used to store the parameters desired by the user. It requires a JSON file to be able to funciton, as it reads directly from the file and stores the values in a variable.
    
    """

    def __init__(self, json_file):

        with open(json_file, 'r') as json_file:
            read_param = json.load(json_file)

        if 'Wire' in read_param:
            self.type_S = 'wire'
            self.B = read_param['Wire']['B']
            self.C1 = read_param['Wire']['C1']
            self.N = read_param['Wire']['N']
            self.M = read_param['Wire']['M']
            self.length = read_param['Wire']['Length']
            self.radius = read_param['Wire']['Radius']
            self.ns = read_param['Wire']['Number_segments']
            self.height = 0
            self.width = 0
            self.raw_stl = read_param['Wire']['Raw_stl']
            self.nfaces = 0
            self.surfaces = 0
            self.energies = 0
            self.n_at = 0


        elif 'Box' in read_param:
            self.type_S = 'box'
            self.N = read_param['Box']['N']
            self.M = read_param['Box']['M']
            self.B = read_param['Box']['B']
            self.C1 = read_param['Box']['C1']
            self.height = read_param['Box']['Height']
            self.length = read_param['Box']['Length']
            self.width = read_param['Box']['Width']
            self.ns = read_param['Box']['Number_segments']
            self.radius = 0
            self.raw_stl = read_param['Box']['Raw_stl']
            self.nfaces = 0
            self.surfaces = 0
            self.energies = 0
            self.n_at = 0


        elif 'Sphere' in read_param:
            self.type_S = 'sphere'
            self.B = read_param['Sphere']['B']
            self.N = read_param['Sphere']['N']
            self.M = read_param['Sphere']['M']
            self.C1 = read_param['Sphere']['C1']
            self.radius = read_param['Sphere']['Radius']
            self.ns = read_param['Sphere']['Number_segments']
            self.height = 0
            self.length = 0
            self.width = 0
            self.raw_stl = read_param['Sphere']['Raw_stl']
            self.nfaces = 0
            self.surfaces = 0
            self.energies = 0
            self.n_at = 0

        elif 'Poly' in read_param:
            self.type_S = 'poly'
            self.B = read_param['Poly']['B']
            self.C1 = read_param['Poly']['C1']
            self.N = read_param['Poly']['N']
            self.M = read_param['Poly']['M']
            self.length = read_param['Poly']['Length']
            self.nfaces = read_param['Poly']['N_Faces']
            self.radius = read_param['Poly']['Radius']
            self.ns = read_param['Poly']['Number_segments']
            self.height = 0
            self.width = 0
            self.raw_stl = read_param['Poly']['Raw_stl']
            self.surfaces = 0
            self.energies = 0
            self.n_at = 0

        elif 'Wulff' in read_param:
            self.type_S = 'wulff'
            self.B = read_param['Wulff']['B']
            self.C1 = read_param['Wulff']['C1']
            self.n_at = read_param['Wulff']['n_atoms']
            self.surfaces = read_param['Wulff']['Surfaces']
            self.energies = read_param['Wulff']['Energies']
            self.N = read_param['Wulff']['N']
            self.M = read_param['Wulff']['M']
            self.length = 0
            self.nfaces = 0
            self.radius = 0
            self.ns = 0
            self.height = 0
            self.width = 0
            self.raw_stl = read_param['Wulff']['Raw_stl']

        elif 'Cube' in read_param:
            self.type_S = 'cube'
            self.B = read_param['Cube']['B']
            self.C1 = read_param['Cube']['C1']
            self.n_at = 0
            self.surfaces = 0
            self.energies = 0
            self.N = read_param['Cube']['N']
            self.M = read_param['Cube']['M']
            self.length = read_param['Cube']['Length']
            self.nfaces = 0
            self.radius = 0
            self.ns = 0
            self.height = read_param['Cube']['Length']
            self.width = read_param['Cube']['Length']
            self.raw_stl = read_param['Cube']['Raw_stl']

        # This loop is storing all the parameters from the json file and storing them, so they can be called in the
        # main code to fun the functions
        self.lattice_structure = read_param['ATOM_Param']['Lattice_structure']
        self.lattice_parameter = read_param['ATOM_Param']['Lattice_parameter']
        self.material = read_param['ATOM_Param']['Material']
        self.orien_x = read_param['ATOM_Param']['Orien_x']
        self.orien_y = read_param['ATOM_Param']['Orien_y']
        self.orien_z = read_param['ATOM_Param']['Orien_z']

    def output(self, json_file):
        """
        This function obtains the desired file type of the user. Having the option to obtain an stl file or lmp file based on the decsion if MD or FEM is desired.

        :param json_file: A json file that that has the parameters stored into it
        :type json_file: file
        """
        with open(json_file, 'r') as f:
            read_param = json.load(f)
        if read_param['Output']['FEM_stl'] == 1:
            return ('FEM_stl')
        # creates a mesh
        # can call from the other functions and get the value of the parameters
        elif read_param['Output']['ATOM_lmp'] == 1:
            return ('ATOM_lmp')
        # create a mesh
        # then that file continues in the code
        elif read_param['Output']['ATOM_lmp'] == 0:
            pass
        elif read_param['Output']['FEM_stl'] == 0:
            pass
