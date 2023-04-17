# ---------------------------------------------------------------------------
# 
# Title: SAMPLE CLASS
#
# Authors: Jonathan Amodeo, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot & Hugo Iteney
#
# Date: June 01, 2022
#
# The SAMPLE_class is a class that obtains the desired object shape (nanowire,slab,sphere) 
# and creates a desired file type. MAKE_STL creates an stl file of the shape, this
# stl file is used for FEM. In this def the creation of the shape is based on user 
# desired parameters that were acquired through the PARAM_class.
# The MAKE_ATOM creates a lmp file that is used for atomsk. It requires the stl made in the MAKE_STL def.
# The parameters for the creation of a lmp file were based of the desired parameters requested by the user.
# ---------------------------------------------------------------------------

from src import Func_pyrough as fp
import numpy as np
import subprocess

# np.set_printoptions(threshold=sys.maxsize)

class Sample(object):
    def __init__(self, type_sample):
        self.sample = type_sample

    def make_stl(self,
                 type_sample,
                 H,
                 C1,
                 RMS,
                 N,
                 M,
                 radius,
                 length,
                 height,
                 width,
                 ns,
                 raw_stl,
                 nfaces,
                 surfaces,
                 energies,
                 n_at,
                 lattice_structure,
                 lattice_parameter,
                 material,
                 orien_x,
                 orien_y,
                 orien_z,
                 out_pre
                 ):
        """
        Creates an stl file based on the parameters chosen by the user.
        Based on the type of sample entered, the code will sort through the options of type of samples and if true then the code under the type of sample will be excuted.
        The type of sample affects the type of stl file returned.
        The code executed requires the parameters to create the object stl and to apply surface roughness onto the object. ...;
        In continuation will add that a separate files will be stored if they want to redo the code.

        :param type_sample: Type of the sample
        :type type_sample: str
        :param B: The degree the roughness is dependent on
        :type B: float
        :param C1: Roughness normalization factor
        :type C1: float
        :param N: Length of the random numbers matrix
        :type N: int
        :param M: Length of the random numbers matrix
        :type M: int
        :param radius: The radius of the nanowire or the sphere
        :type radius: float
        :param length: The length of the box
        :type length: float
        :param height: The height of the box
        :type height: float
        :param width: The width of the box
        :type width: float
        :param raw_stl: Input STL file
        :type raw_stl: stl
        :param nfaces: Number of faces in the case of a faceted wire
        :type nfaces: int
        :param surfaces: List of surfaces used for Wulff theory
        :type surfaces: list
        :param energies: List of energies for each associated surface
        :type energies: list
        :param n_at: Number of atoms in the case of a faceted NP
        :type n_at: int
        :param lattice_structure: Crystallographic structure of the material
        :type lattice_structure: str
        :param lattice_parameter: Lattice parameter of the material
        :type lattice_parameter: float
        :param material: Type of material
        :type material: str

        :returns: List of nodes and stl file name
        """
        if type_sample == 'wire':
            vertices, stl = make_wire(type_sample,
                                      2*(1+H),
                                      C1,
                                      RMS,
                                      N,
                                      M,
                                      radius,
                                      length,
                                      ns,
                                      raw_stl,
                                      out_pre
                                      )

            return (vertices, stl)

        elif type_sample == 'box':
            vertices, stl = make_box(type_sample,
                                     2*(1+H),
                                     C1,
                                     RMS,
                                     N,
                                     M,
                                     length,
                                     height,
                                     width,
                                     ns,
                                     raw_stl,
                                     out_pre
                                     )

            return (vertices, stl)

        elif type_sample == 'sphere':
            vertices, stl = make_sphere(type_sample,
                                        4*H,
                                        C1,
                                        RMS,
                                        N,
                                        radius,
                                        ns,
                                        raw_stl,
                                        out_pre
                                        )

            return (vertices, stl)

        elif type_sample == 'poly':
            vertices, stl = make_poly(type_sample,
                                      2*(1+H),
                                      C1,
                                      RMS,
                                      N,
                                      M,
                                      length,
                                      nfaces,
                                      radius,
                                      ns,
                                      raw_stl,
                                      out_pre,
                                      )
            return (vertices, stl)

        elif type_sample == 'wulff':
            vertices, stl = make_wulff(type_sample,
                                       2*(1+H),
                                       C1,
                                       RMS,
                                       N,
                                       M,
                                       surfaces,
                                       energies,
                                       n_at,
                                       ns,
                                       raw_stl,
                                       lattice_structure,
                                       lattice_parameter,
                                       material,
                                       orien_x,
                                       orien_y,
                                       orien_z,
                                       out_pre
                                       )

            return (vertices, stl)

        elif type_sample == 'cube':
            vertices, stl = make_cube(type_sample,
                                      2*(1+H),
                                      C1,
                                      RMS,
                                      N,
                                      M,
                                      length,
                                      ns,
                                      raw_stl,
                                      out_pre
                                      )

            return (vertices, stl)

    def make_atom(self,
                  STL,
                  lattice_str,
                  lattice_par,
                  material,
                  orien_x,
                  orien_y,
                  orien_z,
                  vertices,
                  out_pre
                  ):

        """
        Creates sample_with_atoms.lmp which is an atomic position file. This requires the stl of the object that already has surface roughness applied to it. It adds then atoms in it.

        :param STL: The stl file with an object that has surface roughness applied to it
        :type STL: str
        :param lattice_str: The lattice structure of the cristal
        :type lattice_str: str
        :param lattice_par: The lattice parameter of the cristal in Angstrom
        :type lattice_par: float
        :param material: The chemical symbol of the desired element
        :type material: str
        :param orien_x: The orientation of the cristal in the x direction
        :type orien_x: list
        :param orien_y: The orientation of the cristal in the y direction
        :type orien_y: list
        :param orien_z: The orientation of the cristal in the z direction
        :type orien_z: list
        :param vertices: List of nodes
        :type vertices: array
        """
        dim_x = max(vertices[:, 0]) - min(vertices[:, 0])
        dim_y = max(vertices[:, 1]) - min(vertices[:, 1])
        dim_z = max(vertices[:, 2]) - min(vertices[:, 2])
        dup_x, orien_x = fp.duplicate(dim_x, orien_x, lattice_par)
        dup_y, orien_y = fp.duplicate(dim_y, orien_y, lattice_par)
        dup_z, orien_z = fp.duplicate(dim_z, orien_z, lattice_par)

        lattice_par = str(lattice_par)

        subprocess.call(['atomsk', '--create', lattice_str, lattice_par, material, 'orient', orien_x, orien_y, orien_z,
                         '-duplicate', dup_x, dup_y, dup_z, 'material_supercell.lmp'])
        subprocess.call(
            ['atomsk', 'material_supercell.lmp', '-select', 'stl', 'center', STL, '-select', 'invert', '-rmatom',
             'select', out_pre + '.lmp'])
        subprocess.call(['rm', 'material_supercell.lmp'])
        fp.rebox(out_pre + '.lmp')


def make_wire(type_sample,
              B,
              C1,
              RMS,
              N,
              M,
              radius,
              length,
              ns,
              raw_stl,
              out_pre
              ):
    """
    Creates an stl file of a rough wire

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param radius: Radius of the wire
    :type radius: float
    :param length: Length of the box
    :type length: float
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str

    :return: List of nodes and STL file name
    """
    vertices, faces = fp.read_stl(type_sample, raw_stl, 0, length, 0, radius, ns, 0,
                                  out_pre)  # reads if the user has inputted a stl file or if the mesh needs to me generated and returnes the faces and the vertices of the stl file

    vertices, nodenumber = fp.node_indexing(
        vertices)  # creates a column that has an assigned index number for each row in vertices; returns vertices with the addtion of this new column and also returns the nodenumbers which is an array file from 0 - length of the vertices

    nodesurf = fp.node_surface(type_sample, vertices, nodenumber, 0,
                               0)  # nodes at the surface of the object. These nodes will have the surface roughness applied to it.

    cy_nodesurf = np.array(
        [fp.rho(nodesurf[:, 0], nodesurf[:, 1]), fp.theta(nodesurf[:, 0], nodesurf[:, 1]), nodesurf[:, 2],
         nodesurf[:, 3]]).T  # CONVERT TO CYLINDRICS COORDINATES (Rho Theta Z nodenumb)

    cy_nodesurf = cy_nodesurf[np.lexsort((cy_nodesurf[:, 1], cy_nodesurf[:, 2]))]

    # creating X and Y vector  (Connection matrix size) going to make one line code
    xv = (cy_nodesurf[:, 1] / max(cy_nodesurf[:, 1]) + 1) / 2
    yv = (cy_nodesurf[:, 2] / max(cy_nodesurf[:, 2]) + 1) / 2

    sfrN, sfrM = fp.vectors(N, M)  # creating vectors for M and N
    m, n = fp.random_numbers(sfrN, sfrM)  # normal gaussian for the amplitude

    z = fp.random_surf2(type_sample, m, n, N, M, B, xv, yv, sfrM, sfrN, C1, RMS,
                        out_pre)  # Returns an array with the Z values that will replace the previous z vlaues in the vertices array these represent the rougness on the surface

    vertices = fp.make_rough(type_sample, z, cy_nodesurf, vertices, 0)

    vertices = vertices[:,
               :3]  # gets ride of the index column because stl file generator takes only a matrix with 3 columns

    fp.stl_file(vertices, faces, out_pre)  # creates an stl file of the cylinder with roughness on the surface
    return (vertices, out_pre + '.stl')  # returns the stl file name


def make_box(type_sample,
             B,
             C1,
             RMS,
             N,
             M,
             length,
             height,
             width,
             ns,
             raw_stl,
             out_pre
             ):
    """
    Creates an stl file of a rough box

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param length: Length of the box
    :type length: float
    :param height: Height of the box
    :type height: float
    :param width: Width of the box
    :type width: float
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str

    :return: List of nodes and STL file name
    """
    vertices, faces = fp.read_stl(type_sample, raw_stl, width, length, height, 0, ns, 0,
                                  out_pre)  # reads if the user has inputted a stl file or if the mesh needs to me generated and returnes the faces and the vertices of the stl file
    vertices, nodenumber = fp.node_indexing(
        vertices)  # creates a column that has an assigned index number for each row in vertices; returns vertices with the addtion of this new column and also returns the nodenumbers which is an array fille from 0 - length of the vertices

    nodesurf = fp.node_surface(type_sample, vertices, nodenumber, 0,
                               0)  # nodes at the surface of the object. These nodes will have the surface roughness applied to it.
    nodesurf = nodesurf[np.lexsort((nodesurf[:, 0], nodesurf[:, 1]))]

    xv = nodesurf[:, 0] / max(nodesurf[:, 0])
    yv = nodesurf[:, 1] / max(nodesurf[:, 1])

    sfrN, sfrM = fp.vectors(N, M)  # creating vectors for M and N
    m, n = fp.random_numbers(sfrN, sfrM)  # normal gaussian for the amplitude

    z = fp.random_surf2(type_sample, m, n, N, M, B, xv, yv, sfrM, sfrN, C1, RMS, out_pre)  # Returns an array with the Z values that will replace the previous z vlaues in the vertices array these represent the rougness on the surface
    vertices = fp.make_rough(type_sample, z, nodesurf, vertices, 0)

    vertices = vertices[:,
               :3]  # gets rid of the index column because stl file generator takes only a matrix with 3 columns

    fp.stl_file(vertices, faces, out_pre)  # creates an stl file of the box with roughness on the surface
    return (vertices, out_pre + '.stl')  # returns the stl file name


def make_sphere(type_sample,
                B,
                C1,
                RMS,
                N,
                radius,
                ns,
                raw_stl,
                out_pre
                ):
    """
    Creates an stl file of a rough sphere

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param radius: Radius of the sphere
    :type radius: float
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str

    :return: List of nodes and STL file name
    """
    vertices, faces = fp.read_stl(type_sample, raw_stl, 0, 0, 0, radius, ns, 0,
                                  out_pre)  # reads if the user has inputted a stl file or if the mesh needs to me generated and returnes the faces and the vertices of the stl file

    nbPoint = len(vertices)  # Stores an int that repesents the number of points on the sphere

    r = np.zeros(
        nbPoint)  # an array with only zeros, with the amount of zeros being equal to the number of points in the sphere

    x, y, z, t = fp.coord_sphere(
        vertices)  # Creats a matirx with the columns that correspond to the coordinates of either x, y, z and t which creates a matrix of the type

    vert_phi_theta = fp.vertex_tp(x, y, t,
                                  z)  # Creates an array filled with two elements that are the angles coresponding to the postion of the node on the sphere.

    thetaa = fp.theta(y, x)  # Calculates the arctan2 of two given arrays whose size are the same.
    phii = fp.phi(t, z)  # Calculates the arctan2 of an array filled with vector norms and an arrray filled with z cooridnates which are the same size.

    r = fp.rough_matrix_sphere(nbPoint, B, thetaa, phii, vert_phi_theta, C1, RMS, r)  # creates the displacement values of the nodes on the surface of the sphere

    C2 = 1. / nbPoint / N / 2.

    fp.stat_sphere(r, C1, C2, out_pre)  # prints the statistics of the sphere

    new_vertex = fp.coord_cart_sphere(C1, C2, r, vertices, t, z, y, x)  # creates a new matrix with x, y, z in cartesian coordinates

    fp.stl_file(new_vertex, faces, out_pre)  # creates an stl file of the sphere with roughness on the surface

    return (vertices, out_pre + '.stl')  # returns the stl file name


def make_poly(type_sample,
              B,
              C1,
              RMS,
              N,
              M,
              length,
              nfaces,
              radius,
              ns,
              raw_stl,
              out_pre,
              ):
    """
    Creates an stl file of a faceted wire

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param length: Length of the wire
    :type length: float
    :param nfaces: Number of faces of the wire
    :type nfaces: int
    :param radius: Radius of the wire
    :type radius: float
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str

    :return: List of nodes and STL file name
    """
    points, angles = fp.base(radius, nfaces)

    vertices, faces = fp.read_stl(type_sample, raw_stl, 0, length, 0, radius, ns, points,
                                  out_pre)  # reads if the user has inputted a stl file or if the mesh needs to me generated and returnes the faces and the vertices of the stl file

    vertices, nodenumber = fp.node_indexing(
        vertices)  # creates a column that has an assigned index number for each row in vertices; returns vertices with the addtion of this new column and also returns the nodenumbers which is an array fille from 0 - length of the vertices

    nodesurf = fp.node_surface(type_sample, vertices, nodenumber, points,
                               0)  # nodes at the surface of the object. These nodes will have the surface roughness applied to it.

    cy_nodesurf = fp.cart2cyl(nodesurf)
    sorted_cy_nodesurf = cy_nodesurf[np.lexsort((cy_nodesurf[:, 1], cy_nodesurf[:, 2]))]
    nodesurf = fp.cyl2cart(sorted_cy_nodesurf)

    Z = len(np.where(sorted_cy_nodesurf[:, 2] == 0)[0])
    T = len(np.unique(sorted_cy_nodesurf[:, 2]))
    xv = np.linspace(0, 1, Z)
    yv = np.linspace(0, 1, T)
    xv, yv = np.meshgrid(xv, yv)

    sfrN, sfrM = fp.vectors(N, M)  # creating vectors for M and N
    m, n = fp.random_numbers(sfrN, sfrM)  # normal gaussian for the amplitude

    z = fp.random_surf2(type_sample, m, n, N, M, B, xv, yv, sfrM, sfrN, C1, RMS,
                        out_pre)  # Returns an array with the Z values that will replace the previous z vlaues in the vertices array these represent the rougness on the surface

    vertices = fp.make_rough(type_sample, z, nodesurf, vertices, angles)

    vertices = vertices[:,
               :3]  # gets rid of the index column because stl file generator takes only a matrix with 3 columns
    vertices = fp.align_poly(vertices, angles)
    fp.stl_file(vertices, faces, out_pre)  # creates an stl file of the box with roughness on the surface

    return (vertices, out_pre + '.stl')


def make_wulff(type_sample,
               B,
               C1,
               RMS,
               N,
               M,
               surfaces,
               energies,
               n_at,
               ns,
               raw_stl,
               lattice_structure,
               lattice_parameter,
               material,
               orien_x,
               orien_y,
               orien_z,
               out_pre
               ):
    """
    Creates an stl file of a Wulff-shaped NP

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param surfaces: Orientation of surfaces used for Wulff theory
    :type surfaces: list
    :param energies: Energies of associated surfaces for Wulff theory
    :type energies: list
    :param n_at: Number of atoms
    :type n_at: int
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str
    :param lattice_structure: The lattice structure of the cristal
    :type lattice_structure: str
    :param lattice_parameter: The lattice parameter of the cristal in Angstrom
    :type lattice_parameter: float
    :param material: The chemical symbol of the desired element
    :type material: str
    :param orien_x: Orientation along x-axis
    :type orien_x: list
    :param orien_y: Orientation along y-axis
    :type orien_y: list
    :param orien_z: Orientation along z-axis
    :type orien_z: list
    :param out_pre: Prefix for output files
    :type out_pre: str

    :return: List of nodes and STL file name
    """
    obj_points, obj_faces = fp.make_obj(surfaces, energies, n_at, lattice_structure, lattice_parameter, material, orien_x, orien_y, orien_z, out_pre)
    subprocess.call(['rm', out_pre + '.obj'])
    vertices, faces = fp.read_stl_wulff(raw_stl, obj_points, obj_faces,
                                        out_pre)  # reads if the user has inputted a stl file or if the mesh needs to me generated and returnes the faces and the vertices of the stl file
    list_n = fp.faces_normals(obj_points, obj_faces)

    vertices, nodenumber = fp.node_indexing(
        vertices)  # creates a column that has an assigned index number for each row in vertices; returns vertices with the addtion of this new column and also returns the nodenumbers which is an array fille from 0 - length of the vertices

    nodesurf = fp.node_surface(type_sample, vertices, nodenumber, obj_points, obj_faces)
    node_edge, node_corner = fp.node_corner(nodesurf)

    vertices = fp.make_rough_wulff(vertices, B, C1, RMS, N, M, nodesurf, node_edge, node_corner, list_n)

    vertices = vertices[:, :3]

    fp.stl_file(vertices, faces, out_pre)

    return (vertices, out_pre + '.stl')


def make_cube(type_sample,
              B,
              C1,
              RMS,
              N,
              M,
              length,
              ns,
              raw_stl,
              out_pre
              ):
    """
    Creates an stl file of a cubic NP

    :param type_sample: The name of the sample
    :type type_sample: str
    :param B: The degree the roughness is dependent on
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param length: Dimension of the cube
    :type length: float
    :param ns: The number of segments desired
    :type ns: int
    :param raw_stl: Name of the input stl file
    :type raw_stl: str

    :return: List of nodes and STL file name
    """

    obj_points, obj_faces = fp.cube_faces(length)
    list_n = fp.faces_normals(obj_points, obj_faces)

    vertices, faces = fp.read_stl(type_sample, raw_stl, 0, length, 0, 0, ns, 0, out_pre)
    vertices, nodenumber = fp.node_indexing(vertices)

    nodesurf = fp.node_surface(type_sample, vertices, nodenumber, obj_points, obj_faces)
    node_edge, node_corner = fp.node_corner(nodesurf)

    vertices = fp.make_rough_wulff(vertices, B, C1, RMS, N, M, nodesurf, node_edge, node_corner, list_n)

    vertices = vertices[:, :3]

    fp.stl_file(vertices, faces, out_pre)

    return (vertices, out_pre + '.stl')

def make_atom_grain(STL,
                    lattice_structure1,
                    lattice_parameter1,
                    material1,
                    orien_x1,
                    orien_y1,
                    orien_z1,
                    lattice_structure2,
                    lattice_parameter2,
                    material2,
                    orien_x2,
                    orien_y2,
                    orien_z2,
                    vertices,
                    out_pre):
    """
    Creates grain.lmp which is an atomic position file.
    This requires the stl of the box object that already has surface roughness applied to it.
    It adds then atoms of material1 in it and fills the rest of the box with atoms of material2.

    :param STL: The stl file with an object that has surface roughness applied to it
    :type STL: str
    :param lattice_structure1: The lattice structure of the cristal 1
    :type lattice_structure1: str
    :param lattice_parameter1: The lattice parameter of the cristal 1 in Angstrom
    :type lattice_parameter1: float
    :param material1: The chemical symbol of the desired element 1
    :type material1: str
    :param orien_x1: The orientation of the cristal 1 in the x direction
    :type orien_x1: list
    :param orien_y1: The orientation of the cristal 1 in the y direction
    :type orien_y1: list
    :param orien_z1: The orientation of the cristal 1 in the z direction
    :type orien_z1: list
    :param lattice_structure2: The lattice structure of the cristal 2
    :type lattice_structure2: str
    :param lattice_parameter2: The lattice parameter of the cristal 2 in Angstrom
    :type lattice_parameter2: float
    :param material2: The chemical symbol of the desired element 2
    :type material2: str
    :param orien_x2: The orientation of the cristal 2 in the x direction
    :type orien_x2: list
    :param orien_y2: The orientation of the cristal 2 in the y direction
    :type orien_y2: list
    :param orien_z2: The orientation of the cristal 2 in the z direction
    :type orien_z2: list
    :param vertices: List of nodes
    :type vertices: array
    :param out_pre: prefix for output file
    :type out_pre: str

    :return:
    """

    dim_x = max(vertices[:, 0]) - min(vertices[:, 0])
    dim_y = max(vertices[:, 1]) - min(vertices[:, 1])
    dim_z = max(vertices[:, 2]) - min(vertices[:, 2])
    dup_x1, orien_x1 = fp.duplicate(dim_x, orien_x1, lattice_parameter1)
    dup_y1, orien_y1 = fp.duplicate(dim_y, orien_y1, lattice_parameter1)
    dup_z1, orien_z1 = fp.duplicate(2*dim_z, orien_z1, lattice_parameter1)
    dup_x2, orien_x2 = fp.duplicate(dim_x, orien_x2, lattice_parameter2)
    dup_y2, orien_y2 = fp.duplicate(dim_y, orien_y2, lattice_parameter2)
    dup_z2, orien_z2 = fp.duplicate(2*dim_z, orien_z2, lattice_parameter2)

    lattice_parameter1 = str(lattice_parameter1)
    lattice_parameter2 = str(lattice_parameter2)

    subprocess.call(['atomsk', '--create', lattice_structure1, lattice_parameter1, material1, 'orient', orien_x1, orien_y1, orien_z1,
                     '-duplicate', dup_x1, dup_y1, dup_z1, 'mat1_supercell.atsk'])
    subprocess.call(['atomsk', 'mat1_supercell.atsk', '-select', 'stl', STL, '-rmatom', 'select', 'mat1_out.atsk'])

    subprocess.call(['atomsk', '--create', lattice_structure2, lattice_parameter2, material2, 'orient', orien_x2, orien_y2, orien_z2,
                     '-duplicate', dup_x2, dup_y2, dup_z2, 'mat2_supercell.atsk'])
    subprocess.call(['atomsk', 'mat2_supercell.atsk', '-select', 'stl', STL, '-select', 'invert', '-rmatom', 'select', 'mat2_out.atsk'])

    subprocess.call(['atomsk', '--merge', '2', 'mat1_out.atsk', 'mat2_out.atsk', out_pre+'.lmp'])

    subprocess.call(['rm', 'mat1_supercell.atsk', 'mat2_supercell.atsk', 'mat1_out.atsk', 'mat2_out.atsk'])
    #fp.rebox(out_pre + '.lmp')
