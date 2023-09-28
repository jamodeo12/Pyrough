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

        if 'Grain' in read_param :
            self.type_S = 'grain'
            self.N = read_param['Grain']['N']
            self.M = read_param['Grain']['M']
            self.H = read_param['Grain']['H']
            try:
                self.C1 = read_param['Grain']['C1']
            except KeyError:
                self.C1 = ''
            try:
                self.RMS = read_param['Grain']['RMS']
            except KeyError:
                self.RMS = ''
            self.height = read_param['Grain']['Height']
            self.length = read_param['Grain']['Length']
            self.width = read_param['Grain']['Width']
            self.ns = read_param['Grain']['Mesh size']
            self.raw_stl = read_param['Grain']['Raw_stl']
            self.lattice_structure1 = read_param['ATOM1_Param']['Lattice_structure']
            self.lattice_parameter1 = read_param['ATOM1_Param']['Lattice_parameter']
            self.material1 = read_param['ATOM1_Param']['Material']
            self.orien_x1 = read_param['ATOM1_Param']['Orien_x']
            self.orien_y1 = read_param['ATOM1_Param']['Orien_y']
            self.orien_z1 = read_param['ATOM1_Param']['Orien_z']
            self.lattice_structure2 = read_param['ATOM2_Param']['Lattice_structure']
            self.lattice_parameter2 = read_param['ATOM2_Param']['Lattice_parameter']
            self.material2 = read_param['ATOM2_Param']['Material']
            self.orien_x2 = read_param['ATOM2_Param']['Orien_x']
            self.orien_y2 = read_param['ATOM2_Param']['Orien_y']
            self.orien_z2 = read_param['ATOM2_Param']['Orien_z']

        else :
            if 'cWire' in read_param:
                self.type_S = 'wire'
                self.H = read_param['cWire']['H']
                try:
                    self.C1 = read_param['cWire']['C1']
                except KeyError:
                    self.C1 = ''
                try:
                    self.RMS = read_param['cWire']['RMS']
                except KeyError:
                    self.RMS = ''
                self.N = read_param['cWire']['N']
                self.M = read_param['cWire']['M']
                self.length = read_param['cWire']['Length']
                self.radius = read_param['cWire']['Radius']
                self.ns = read_param['cWire']['Mesh_size']
                try:
                    self.alpha = read_param['cWire']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param['cWire']['Raw_stl']
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']


            elif 'Box' in read_param:
                self.type_S = 'box'
                self.N = read_param['Box']['N']
                self.M = read_param['Box']['M']
                self.H = read_param['Box']['H']
                try:
                    self.C1 = read_param['Box']['C1']
                except KeyError:
                    self.C1 = ''
                try:
                    self.RMS = read_param['Box']['RMS']
                except KeyError :
                    self.RMS = ''
                self.height = read_param['Box']['Height']
                self.length = read_param['Box']['Length']
                self.width = read_param['Box']['Width']
                self.ns = read_param['Box']['Mesh_size']
                try:
                    self.alpha = read_param['Box']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.radius = 0
                self.raw_stl = read_param['Box']['Raw_stl']
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']


            elif 'Sphere' in read_param:
                self.type_S = 'sphere'
                self.H = read_param['Sphere']['H']
                self.N = read_param['Sphere']['N']
                self.M = read_param['Sphere']['M']
                self.C1 = read_param['Sphere']['C1']
                self.RMS = ''
                self.radius = read_param['Sphere']['Radius']
                self.ns = read_param['Sphere']['Mesh_size']
                try:
                    self.alpha = read_param['Sphere']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.length = 0
                self.width = 0
                self.raw_stl = read_param['Sphere']['Raw_stl']
                self.nfaces = 0
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']

            elif 'fWire' in read_param:
                self.type_S = 'poly'
                self.H = read_param['fWire']['H']
                try:
                    self.C1 = read_param['fWire']['C1']
                except KeyError:
                    self.C1 = ''
                try:
                    self.RMS = read_param['fWire']['RMS']
                except KeyError:
                    self.RMS = ''
                self.N = read_param['fWire']['N']
                self.M = read_param['fWire']['M']
                self.length = read_param['fWire']['Length']
                self.nfaces = read_param['fWire']['N_Faces']
                self.radius = read_param['fWire']['Radius']
                self.ns = read_param['fWire']['Mesh_size']
                try:
                    self.alpha = read_param['fWire']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param['fWire']['Raw_stl']
                self.surfaces = 0
                self.energies = 0
                self.n_at = 0
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']

            elif 'Wulff' in read_param:
                self.type_S = 'wulff'
                self.H = read_param['Wulff']['H']
                try:
                    self.C1 = read_param['Wulff']['C1']
                except KeyError:
                    self.C1 = ''
                try:
                    self.RMS = read_param['Wulff']['RMS']
                except KeyError:
                    self.RMS = ''
                self.n_at = read_param['Wulff']['n_atoms']
                self.surfaces = read_param['Wulff']['Surfaces']
                self.energies = read_param['Wulff']['Energies']
                self.N = read_param['Wulff']['N']
                self.M = read_param['Wulff']['M']
                self.length = 0
                self.nfaces = 0
                self.radius = 0
                self.ns = read_param['Wulff']['Mesh_size']
                try:
                    self.alpha = read_param['Wulff']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.height = 0
                self.width = 0
                self.raw_stl = read_param['Wulff']['Raw_stl']
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']

            elif 'Cube' in read_param:
                self.type_S = 'cube'
                self.H = read_param['Cube']['H']
                try:
                    self.C1 = read_param['Cube']['C1']
                except KeyError:
                    self.C1 = ''
                try:
                    self.RMS = read_param['Cube']['RMS']
                except KeyError:
                    self.RMS = ''
                self.n_at = 0
                self.surfaces = 0
                self.energies = 0
                self.N = read_param['Cube']['N']
                self.M = read_param['Cube']['M']
                self.length = read_param['Cube']['Length']
                self.nfaces = 0
                self.radius = 0
                self.ns = read_param['Cube']['Mesh_size']
                try:
                    self.alpha = read_param['Cube']['Refine_factor']
                except KeyError:
                    self.alpha = 1
                self.height = read_param['Cube']['Length']
                self.width = read_param['Cube']['Length']
                self.raw_stl = read_param['Cube']['Raw_stl']
                self.ext_fem = read_param['Output']['FEM']
                self.ext_ato = read_param['Output']['ATOM']

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
        if len(read_param['Output']['ATOM']) > 0:
            return ('ATOM')

