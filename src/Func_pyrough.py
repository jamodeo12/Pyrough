# ---------------------------------------------------------------------------
# 
# Title: Func_pyrough
#
# Authors: Jonathan Amodeo, Javier Gonzalez, Jennifer Izaguirre, Christophe Le Bourlot & Hugo Iteney
#
# Date: June 01, 2022
#
# 
# Func_pyrough are all the functions used in pyrough. These are required to excute the Pyroughs main code.
#  
# 
# ---------------------------------------------------------------------------

import pygmsh
import meshio
import numpy as np
import sys
import math
import scipy.special as sp
from wulffpack import SingleCrystal
from ase.build import bulk
from pathlib import Path

np.set_printoptions(threshold=sys.maxsize)


def rdnsurf(m, n, B, xv, yv, sfrM, sfrN):
    """
    Generates random surface roughness that will replace the previous z values in the vertices matrix.

    :type m: array 
    :param m: Wavenumber
    :type n: array
    :param n: Wavenumber
    :type B: float
    :param B: The degree the roughness is dependent on 
    :type xv: array
    :param xv: Unique x coordinate values from the objects vertices 
    :type yv: array
    :param yv: Unique y coordinate values from the objects vertices 
    :type sfrM: array
    :param sfrM: Matrix of random numbers that range from 1 to 2N 
    :type sfrN: array
    :param sfrN: Matrix of random numbers that range from 1 to 2N

    :return: Roughness height matrix
    """
    print('====== > Creating random surface....')
    Z = 0.0
    for i in range(len(sfrM)):
        for j in range(len(sfrN)):
            if sfrM[i] == 0 and sfrN[j] == 0:
                continue
            else:
                mod = (sfrM[i] ** 2 + sfrN[j] ** 2) ** (-0.5 * B)
                Z = Z + m[i][j] * mod * np.cos(2 * np.pi * (sfrM[i] * xv + sfrN[j] * yv) + n[i][j])
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter3D(xv, yv, Z, 'gray')
    # # ax.contour3D(xv, yv, Z, cmap=cm.coolwarm, linewidth=0)
    # # ax.plot_surface(xv, yv, Z, rstride=1, cstride=1, cmap=cm.nipy_spectral, linewidth=0, antialiased=False)
    # ax.view_init(elev=90, azim=-90)
    # plt.show()
    return Z


def rho(x, y):
    """
    The Pythagorean theorem equation is used to obtain a value for a side of a triangle. 
    
    :param x: Represents the length of a side in the triangle 
    :type x: int
    :param y: Represents the length of a side in the triangle 
    :type y: int
    
    :return: The length of hypotenuse
    """
    return np.sqrt(x ** 2 + y ** 2)


def theta(x, y):
    """
    Calculates the arctan2 of two given arrays whose size are the same. 

    :param x: x coordinates 
    :type x: array
    :param y: y coordinates
    :type y: array

    :return: An array of angles in radians
    """
    return np.arctan2(y, x)


def cylinder(l, r, ns, out_pre):
    """
    Creates a cylindrical mesh.
    
    :param l: Length of the desired nanowire 
    :type l: float
    :param r: Radius of the nanowire 
    :type r: float
    :param ns: The number of segments desired on the nanowire 
    :type ns: int
    :param out_pre: Prefix for output files
    :type out_pre: str

    :returns: List of points and faces
    """
    print('====== > Creating the Mesh')
    points = []
    theta_list = np.linspace(0, 2 * np.pi, ns)
    for theta in theta_list:
        points.append([r * np.cos(theta), r * np.sin(theta)])
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(points, mesh_size=3)
        geom.extrude(poly, [0, 0, l], num_layers=75)
        mesh = geom.generate_mesh()
    vertices = mesh.points
    faces = mesh.get_cells_type('triangle')
    write_stl('Raw_' + out_pre + '.stl', vertices, faces)
    print('====== > Done creating the Mesh')
    return (vertices, faces)


def box(width, length, height, ns, out_pre):
    """
    Creates a box mesh.

    :param width: Width of the desired box
    :type width: float
    :param length: Length of the desired box
    :type length: float
    :param height: Height of the desired box
    :type height: float
    :param ns: Number of segments
    :type ns: int
    :param out_pre: Prefix for output files
    :type out_pre: str

    :return: List of points and faces
    """
    print('====== > Creating the Mesh')
    with pygmsh.geo.Geometry() as geom:
        geom.add_box(0, length, 0, width, 0, height, mesh_size=4)
        mesh = geom.generate_mesh()
    vertices = mesh.points
    faces = mesh.get_cells_type('triangle')
    mesh.write('Raw_' + out_pre + '.stl')
    print('====== > Done creating the Mesh')
    return (vertices, faces)


def sphere(r, ns, out_pre):
    """
    Creates a sphere mesh.

    :param r: Radius of the desired sphere
    :type r: float
    :param ns: The number of segments desired on the sphere
    :type ns: int
    :param out_pre: Prefix of output files
    :type out_pre: str

    :returns: List of points and faces
    """
    print('====== > Creating the Mesh')
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_ball([0.0, 0.0, 0.0], r, mesh_size=7)
        mesh = geom.generate_mesh()
    vertices = mesh.points
    faces = mesh.get_cells_type('triangle')
    write_stl('Raw_' + out_pre + '.stl', vertices, faces)
    print('====== > Done creating the Mesh')
    return (vertices, faces)


def poly(length, points, out_pre):
    """
    Creates a faceted wire mesh.

    :param length: Length of the wire
    :type length: float
    :param points: Base points of the wire
    :type points: array
    :param out_pre: Prefix of output files
    :type out_pre: str

    :returns: List of points and faces
    """
    print('====== > Creating the Mesh')
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(points, mesh_size=3)
        geom.extrude(poly, [0.0, 0.0, length], num_layers=100)
        mesh = geom.generate_mesh()
    vertices = mesh.points
    faces = mesh.get_cells_type('triangle')
    mesh.write('Raw_' + out_pre + '.stl')
    print('====== > Done creating the Mesh')
    return (vertices, faces)


def wulff(obj_points, obj_faces, out_pre):
    """
    Creates a Wulff-Shaped NP mesh

    :param obj_points: Corners of the NP (from OBJ file)
    :type obj_points: list
    :param obj_faces: Faces of the NP (from OBJ file)
    :type obj_faces: list

    :returns: List of points and faces
    """
    with pygmsh.geo.Geometry() as geom:
        for polygone in obj_faces:
            list_points = []
            for i in polygone:
                list_points.append(obj_points[int(i - 1)])
            list_points = np.asarray(list_points)
            geom.add_polygon(list_points, mesh_size=3)
        mesh = geom.generate_mesh()
        vertices = mesh.points
        faces = mesh.get_cells_type('triangle')
    write_stl('Raw_' + out_pre + '.stl', vertices, faces)
    return (vertices, faces)


def cube(length, out_pre):
    """
    Creates a cube mesh.

    :param length: Length of the desired box
    :type length: float
    :param out_pre: Prefix for output files
    :type out_pre: str

    :return: List of points and faces
    """
    print('====== > Creating the Mesh')
    with pygmsh.geo.Geometry() as geom:
        geom.add_box(0, length, 0, length, 0, length, mesh_size=length / 30)
        mesh = geom.generate_mesh()
    vertices = mesh.points
    faces = mesh.get_cells_type('triangle')
    mesh.write('Raw_' + out_pre + '.stl')
    print('====== > Done creating the Mesh')
    return (vertices, faces)


def make_obj(surfaces, energies, n_at, lattice_structure, lattice_parameter, material, orien_x, orien_y, orien_z,
             out_pre):
    """
    Creates an OBJ file of a faceted NP. Stores the points and faces from this file.

    :param surfaces: List of surfaces for Wulff theory
    :type surfaces: list
    :param energies: List of energies associated to surfaces for Wulff theory
    :type energies: list
    :param n_at: Number of atoms
    :type n_at: int
    :param lattice_structure: Atomic structure of the sample
    :type lattice_structure: str
    :param lattice_parameter: Lattice parameter
    :type lattice_parameter: float
    :param material: Type of atom
    :type material: str
    :param orien_x: Orientation along x-axis
    :type orien_x: list
    :param orien_y: Orientation along y-axis
    :type orien_y: list
    :param orien_z: Orientation along z-axis
    :type orien_z: list
    :param out_pre: Prefix for output files
    :type out_pre: str

    :returns: List of points and faces of OBJ file
    """
    surface_energies = {tuple(surfaces[i]): float(energies[i]) for i in range(0, len(surfaces))}
    prim = bulk(material, crystalstructure=lattice_structure, a=lattice_parameter)
    particle = SingleCrystal(surface_energies, primitive_structure=prim, natoms=n_at)
    particle.write(out_pre + '.obj')
    with open(out_pre + '.obj') as f:
        list_lines = f.readlines()
        obj_points = []
        obj_faces = []
        for line in list_lines:
            if 'v' in line:
                splited = line.split()[1:]
                coord = [float(i) for i in splited]
                obj_points.append(coord)
            if 'f' in line:
                splited = line.split()[1:]
                coord = [float(i) for i in splited]
                obj_faces.append(coord)
    obj_points = np.asarray(obj_points)
    obj_points_f = rotate_obj_wulff(obj_points, orien_x, orien_z)
    return (obj_points_f, obj_faces)


def read_stl(sample_type, raw_stl, width, length, height, radius, ns, points, out_pre):
    """
    Reads an input stl file or creates a new one if no input

    :param sample_type: Name of the sample
    :type sample_type: str
    :param raw_stl: Name of the input stl file
    :type raw_stl: str
    :param width: Width of the box
    :type width: float
    :param length: Length of the box/wire
    :type length: float
    :param height: Height of the box
    :type height: float
    :param radius: Radius of the wire/sphere
    :type radius: float
    :param ns: The number of segments desired
    :type ns: int
    :param points: List of points constituting the base (in case of faceted wire)
    :type points: array
    :param out_pre: Prefix of output files
    :type out_pre: str

    :return: List of points and faces
    """
    vertices = []
    faces = []
    if raw_stl == "na":
        if sample_type == "box" or sample_type == "grain":
            vertices, faces = box(width, length, height, ns, out_pre)
        elif sample_type == "wire":
            vertices, faces = cylinder(length, radius, ns, out_pre)
        elif sample_type == "sphere":
            vertices, faces = sphere(radius, ns, out_pre)
        elif sample_type == "poly":
            vertices, faces = poly(length, points, out_pre)
        elif sample_type == "cube":
            vertices, faces = cube(length, out_pre)
    else:
        mesh = meshio.read(raw_stl)
        vertices, faces = mesh.points, mesh.cells
        faces = faces[0][1]
    return (vertices, faces)


def read_stl_wulff(raw_stl, obj_points, obj_faces, out_pre):
    """
    Reads an input stl file or creates a new one if no input. Wulff case.

    :param raw_stl: Name of the input stl file
    :type raw_stl: str
    :param obj_points: List of points from OBJ file
    :type obj_points: list
    :param obj_faces: List of faces from OBJ file
    :type obj_faces: list
    :param out_pre: Prefix of output files
    :type out_pre: str

    :returns: List of points and faces
    """
    if raw_stl == "na":
        vertices, faces = wulff(obj_points, obj_faces, out_pre)
    else:
        mesh = meshio.read(raw_stl)
        vertices, faces = mesh.vertices, mesh.faces
    return (vertices, faces)


def stl_file(vertices, faces, out_pre):
    """
    Creates an stl file from the vertices and faces of the desired object. 
    
    :param vertices: The coordinates obtained from the mesh 
    :type vertices: array
    :param faces: The faces of the triangles generated from the mesh
    :type vertices: array
    :param out_pre: Prefix of the output files
    :type out_pre: str
    """
    write_stl(out_pre + '.stl', vertices, np.array(faces))
    return


def phi(t, z):
    """
    Calculates the arctan2 of an array filled with vector norms and an array filled with z coordinates which are the same size.
    
    :param t: Coordinates x y in a list
    :type t: array
    :param z: Z coordinates
    :type z: array
    
    :return: An array of angles in radians
    """
    return np.arctan2(np.linalg.norm(t, axis=1), z)


def vertex_tp(x, y, t, z):
    """
    Creates an array filled with two elements that are the angles corresponding to the position of the node on the sphere.
    
    :param x: X coordinates 
    :type x: array
    :param y: Y coordinates 
    :type y: array
    :param t: x & y coordinates stored in a list
    :type t: array
    :param z: Z coordinates
    :type z: array

    :return: An array with angles as elements
    """
    return np.array([theta(y, x), phi(t, z)]).T


def radius(v):
    """
    Calculates the vector norm for the axis 1 of the values in the given array. 

    :param v: The vertices of the sphere 
    :type v: array 
        
    :return: A vector norm
    """
    return np.linalg.norm(v, axis=1)


def stat_analysis(z, N, M, C1, B, sample_type, out_pre):
    """
    Displays the statistical analysis of the surface roughness generator
    
    :param z: Roughness height matrix
    :type z: float
    :param N: Scaling cartesian position 
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param C1: Roughness normalization factor
    :type C1: float
    :param B: The degree the roughness is dependent on
    :type B: float
    :param sample_type: The name of the sample 
    :type sample_type: str
    :param out_pre: Prefix of output files
    :type out_pre: str

    """
    # z_an = np.reshape(z, -1)
    z_an = z.flatten()
    nu_points = len(z_an)
    mean = np.mean(z_an)
    stand = np.std(z_an)
    rms = np.sqrt(np.sum(np.square(z_an)) / len(z_an))
    skewness = np.sum(np.power((z_an - np.mean(z_an)), 3) / len(z_an)) / np.power(np.std(z_an), 3)
    kurtosis = np.sum(np.power((z_an - np.mean(z_an)), 4) / len(z_an)) / np.power(np.std(z_an), 4)

    stats = [N, M, C1, round(0.5 * B - 1, 2), nu_points, mean, stand, rms, skewness, kurtosis]
    stats = list(map(str, stats))
    stats = [sample_type, 'N = ' + stats[0], 'M = ' + stats[1], 'C1 = ' + stats[2], 'H = ' + stats[3],
             'No. points = ' + stats[4],
             'Mean_Value = ' + stats[5], 'Stand_dev = ' + stats[6], 'RMS = ' + stats[7], 'Skewness = ' + stats[8],
             'Kurtosis = ' + stats[9]]

    np.savetxt(out_pre + '_stat.txt', stats, fmt='%s')
    print('')
    print('------------ Random Surface Parameters-----------------')
    print('         N =', N, '  M = ', M, '  C1 = ', C1, '  H = ', round(0.5 * B - 1, 2))
    print('No. points = ', nu_points)
    print('Mean_Value = ', mean)
    print(' Stand_dev = ', stand)
    print('       RMS = ', rms)
    print('  Skewness = ', skewness)
    print('  Kurtosis = ', kurtosis)
    print('--------------------------------------------------------')
    return


def stat_sphere(r, C1, C2, out_pre):
    """
    Displays the statistical analysis of the surface roughness generator with regards to the sphere. 
    
    :param r: Roughness height matrix
    :type r: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param C2: Roughness normalization factor constant for sphere
    :type C2: float
    :param out_pre: Prefix of output files
    :type out_pre: str
    """
    mean = np.mean(r)
    stand = np.std(C1 * C2 * r)
    rms = np.sqrt(C1 * C2 * np.sum(r * r) / len(r))

    stats = [C1, C2, mean, stand, rms]
    stats = list(map(str, stats))
    stats = ['Sphere', 'C1 = ' + stats[0], 'C2 = ' + stats[1], 'Mean_Value = ' + stats[2], 'Stand_dev = ' + stats[3],
             'RMS = ' + stats[4]]

    np.savetxt(out_pre + '_stat.txt', stats, fmt='%s')

    print("ave = {}".format(np.mean(r)))  # out
    print("Stand_dev = {}".format(np.std(C1 * C2 * r)))
    print("RMS = {}".format(np.sqrt(C1 * C2 * np.sum(r * r) / len(r))))
    return


def concatenate_list_data(a_list):
    """
    Combines the elements of a list into one element in the form of a string.
    
    :param a_list: List with elements of type int
    :type a_list: list

    :return: A string
    """
    result = ''
    for element in a_list:
        result += str(element)
    return result


def duplicate(l, orien, lattice_par):
    """
    Takes in a length and an cristal orientation to calculate the duplication factor for atomsk.

    :param l: Length of one of the sides of the object
    :type l: int
    :param orien: Crystal orientaion
    :type orien: list
    :param lattice_par: Lattice parameter
    :type lattice_par: float

    :returns: Duplication factor and string of orientations
    """
    storage = []
    for x in orien:
        squared = x ** 2
        storage.append(squared)
    Total = sum(storage)
    if Total == 6 or Total == 2:
        distance = (lattice_par * (math.sqrt(Total))) / 2
        dup = math.ceil(l / distance)
    else:
        distance = (lattice_par * (math.sqrt(Total)))
        dup = math.ceil(l / distance)
    end_orien = [(concatenate_list_data(orien))]
    x = '[' + "".join([str(i) for i in end_orien]) + ']'
    dup = str(dup)
    return 0.5*distance, dup, x


def random_numbers(sfrN, sfrM):
    """
    Generates the G and U matrices for the mathematical formulation of rough surfaces.

    :param sfrN: Vector for the N decomposition
    :type sfrN: array
    :param sfrM: Vector for the M decomposition
    :type sfrM: array

    :returns: G and U matrices
    """
    m = 0 + 1 * np.random.randn(len(sfrM), len(sfrN))  # Gaussian distributed
    n = -np.pi / 2 + np.pi * np.random.rand(len(sfrM), len(sfrN))  # Uniform distributed
    # np.savetxt('m.txt', m, fmt='%s')
    # np.savetxt('n.txt', n, fmt='%s')
    return m, n


def node_indexing(vertices):
    """
    Creates a column that has an assigned index number for each row in vertices and also the nodenumbers which is an array file from 0 - length of the vertices

    :param vertices: Vector for the points
    :type vertices: array

    :returns: Numbered vertices and raw numbers column.
    """
    nodenumber = range(0, len(vertices))
    vertices = np.insert(vertices, 3, nodenumber, 1)
    return vertices, nodenumber


def node_surface(sample_type, vertices, nodenumber, points, faces):
    """
    Finds the nodes at the surface of the object. These nodes will have the surface roughness applied to it.

    :param sample_type: The name of the sample
    :type sample_type: str
    :param vertices: List of nodes
    :type vertices: array
    :param nodenumber: Number of the corresponding node
    :type nodenumber: array
    :param points : Polygon shape (Faceted wire case)
    :type points: array
    :param faces: Facets list (Wulff case)
    type faces: array

    :return: Surface nodes
    """
    stay = []
    if sample_type == 'wire':
        max_height = max(vertices[:, 1])
        for x in range(0, len(vertices)):
            if rho(vertices[x, 0], vertices[x, 1]) > (max_height - 0.1):
                stay.append(x)
        no_need = np.delete(nodenumber, stay)  # delete from nodenumbers the ones in the surface
        nodesurf = np.delete(vertices, no_need, 0)
        return nodesurf

    elif sample_type == 'box' or sample_type == 'grain':
        eps = 0.000001
        max_height = max(vertices[:, 2])
        for index in range(0, len(vertices)):
            if abs(vertices[index][2] - max_height) <= eps:
                stay.append(index)
        no_need = np.delete(nodenumber, stay)  # delete from nodenumbers the ones in the surface
        nodesurf = np.delete(vertices, no_need, 0)
        return nodesurf

    elif sample_type == 'poly':
        eps = 0.0001
        face = []
        for i in range(len(points)):
            k = i + 1
            if k == len(points):
                k = 0
            p1 = points[i]
            p2 = points[k]
            if p2[0] - p1[0] >= -eps and p2[0] - p1[0] <= eps:
                K = p2[0]
                for verti in vertices:
                    if verti[0] >= K - eps and verti[0] <= K + eps:
                        face = np.append(face, verti, axis=0)
            else:
                A = (p2[1] - p1[1]) / (p2[0] - p1[0])
                B = p2[1] - A * p2[0]
                for verti in vertices:
                    if A * verti[0] + B >= verti[1] - eps and A * verti[0] + B <= verti[1] + eps:
                        face = np.append(face, verti, axis=0)
        l = len(face)
        nodesurf = np.reshape(face, [int(l / 4), 4])
        return nodesurf

    elif sample_type == 'wulff' or sample_type == 'cube':
        nodesurf = []
        for F in faces:
            L = []
            eps = 0.000001
            A = points[int(F[0] - 1)]
            B = points[int(F[1] - 1)]
            C = points[int(F[2] - 1)]
            AB = B - A
            AC = C - A
            for M in vertices:
                AM = M[:3] - A
                matrix = np.array([[AM[0], AB[0], AC[0]], [AM[1], AB[1], AC[1]], [AM[2], AB[2], AC[2]]], dtype=float)
                det = np.linalg.det(matrix)
                if abs(det) <= eps:
                    L.append([M[0], M[1], M[2], M[3]])
            nodesurf.append(L)
        return nodesurf


def vectors(N, M):
    """
    Creates vector of integers between -N and N. Same for M

    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int

    :returns: Vectors
    """
    sfrN = np.linspace(-N, N, 2 * N + 1)
    sfrM = np.linspace(-M, M, 2 * M + 1)
    return sfrN, sfrM


def coord_sphere(vertices):
    """
    Creates a matrix with the columns that correspond to the coordinates of either x, y, z and t which contains x y z coordinates in an array

    :param vertices: List of nodes
    :type vertices: array

    :return: Coordinates
    """
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]
    t = vertices[:, 0:2]
    return (x, y, z, t)


def rough_matrix_sphere(nbPoint, B, thetaa, phii, vert_phi_theta, C1, RMS, r):
    """
    Creates the displacement values of the nodes on the surface of the sphere

    :param nbPoint: Number of points on the sphere
    :type nbPoint: int
    :param B: The degree the roughness is dependent on
    :type B: float
    :param thetaa: Arctan2 of x and y coordinates
    :type thetaa: array
    :param phii: Arctan2 of vector norms and z coordinates
    :type phii: array
    :param vert_phi_theta: Array filled with two elements that are the angles corresponding to the position of the node on the sphere.
    :type vert_phi_theta: array
    :param C1: Roughness normalization factor
    :type C1: float
    :param RMS: Root Mean Square
    :type RMS: float
    :param r: Roughness height matrix
    :type r: int

    :return: Rough matrix
    """
    N_s = 9
    N_e = 17
    for degree in range(N_s, N_e + 1, 1):  # EQUATION
        print("degree: {}".format(degree))
        _r_amplitude = 0 + 1 * np.random.randn(nbPoint)
        _r_phase = -np.pi / 2 + np.pi * np.random.rand(nbPoint)
        mod = degree ** (-B / 2)
        for i, [theta, phi] in enumerate(vert_phi_theta):
            _phase = sp.sph_harm(0, degree, thetaa - theta, phii - phi).real
            _phase = 2 * _phase / _phase.ptp()
            r += _r_amplitude[i] * mod * np.cos(_phase + _r_phase)
    return r


def coord_cart_sphere(C1, C2, r, vertices, t, z, y, x):
    """
    Creates a new matrix with x, y, z in cartesian coordinates

    :param C1: Roughness normalization factor
    :type C1: float
    :param C2: Roughness normalization factor constant for sphere
    :type C2: float
    :param r: Roughness height matrix
    :type r: int
    :param vertices: List of nodes
    :type vertices: array
    :param t: x y z coordinates
    :type t: array
    :param z: z coordinates
    :type z: array
    :param y: y coordinates
    :type y: array
    :param x: x coordinates
    :type x: array

    :return: Cartesian coordinates
    """
    x = (C1 * C2 * r + radius(vertices)) * np.sin(phi(t, z)) * np.cos(theta(x, y))
    y = (C1 * C2 * r + radius(vertices)) * np.sin(phi(t, z)) * np.sin(theta(x, y))
    z = (C1 * C2 * r + radius(vertices)) * np.cos(phi(t, z))
    new_vertex = np.array([x, y, z]).T
    return new_vertex


def rebox(file_lmp):
    """
    Fits the box dimensions with the sample. Also centers the sample.

    :param file_lmp: .lmp file containing the atom positions
    :type file_lmp: str
    """
    fint = open(file_lmp, "r")
    lines = fint.readlines()
    data = np.array([i.split() for i in lines[15:len(lines):1]])
    listN = data[:, 0]
    listi = data[:, 1]
    listx = data[:, 2]
    listy = data[:, 3]
    listz = data[:, 4]
    listN_int = [int(i) for i in listN]
    listi_int = [int(i) for i in listi]
    listx_float = [float(i) for i in listx]
    listy_float = [float(i) for i in listy]
    listz_float = [float(i) for i in listz]
    mx = min(listx_float)
    my = min(listy_float)
    mz = min(listz_float)
    listx_n = [x - mx for x in listx_float]
    listy_n = [y - my for y in listy_float]
    listz_n = [z - mz for z in listz_float]
    compt = 0
    for i in range(0, len(lines), 1):
        if i == 5:
            lines[i] = "{} {} xlo xhi\n".format(math.floor(min(listx_n)), math.ceil(max(listx_n)))
        if i == 6:
            lines[i] = "{} {} ylo yhi\n".format(math.floor(min(listy_n)), math.ceil(max(listy_n)))
        if i == 7:
            lines[i] = "{} {} zlo zhi\n".format(math.floor(min(listz_n)), math.ceil(max(listz_n)))
        if i >= 15:
            lines[i] = "{} {} {} {} {}\n".format(listN_int[compt], listi_int[compt], listx_n[compt], listy_n[compt],
                                                 listz_n[compt])
            compt += 1
    fint.close()
    fend = open(file_lmp, "w")
    fend.writelines(lines)
    fend.close()
    return


def write_stl(filename, vertices, face_list):
    """
    Creates an STL file from faces and vertices.

    :param filename: name of the STL file
    :type filename: str
    :param vertices: list of vertices
    :type vertices: array
    :param face_list: list of faces
    :type face_list: array
    """
    with open(filename, 'w') as f:
        f.write('solid Created by Gmsh \n')
        for face in face_list:
            p1 = vertices[face[0]]
            p2 = vertices[face[1]]
            p3 = vertices[face[2]]
            OA = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]]
            OB = [p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2]]
            normal = np.cross(OA, OB)
            f.write('facet normal {} {} {}\n'.format(normal[0], normal[1], normal[2]))
            f.write('  outer loop\n')
            f.write('    vertex {} {} {}\n'.format(p1[0], p1[1], p1[2]))
            f.write('    vertex {} {} {}\n'.format(p2[0], p2[1], p2[2]))
            f.write('    vertex {} {} {}\n'.format(p3[0], p3[1], p3[2]))
            f.write('  endloop\n')
            f.write('endfacet\n')
        f.write('endsolid Created by Gmsh')
    return


def rms_calc(Z):
    """
    Calculates the RMS of a height distribution

    :param z: height matrix
    :type z: array

    :return: RMS
    """
    Z = np.asarray(Z)
    Z = Z.flatten()
    return (np.sqrt(np.sum(Z * Z) / len(Z)))


def random_surf2(sample_type, m, n, N, M, B, xv, yv, sfrM, sfrN, C1, RMS, out_pre):
    """
    Returns an array with the Z values representing the surface roughness.

    :param sample_type: The type of the sample
    :type sample_type: str
    :param m: Wavenumbers
    :type m: array
    :param n: Wavenumbers
    :type n: array
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param B: The degree of the roughness
    :type B: float
    :param xv: Unique x coordinate values from the objects vertices
    :type xv: array
    :param yv: Unique y coordinate values from the objects vertices
    :type yv: array
    :param sfrN: Vector for the N decomposition
    :type sfrN: array
    :param sfrM: Vector for the M decomposition
    :type sfrM: array
    :param C1: Roughness normalization factor
    :type C1: float
    :param RMS: Root Mean Square
    :type RMS: float
    :param out_pre: Prefix of output files
    :type out_pre: str

    :return: Surface roughness
    """
    Z = rdnsurf(m, n, B, xv, yv, sfrM, sfrN)
    # ax = plt.axes(projection='3d')
    # color_map = plt.get_cmap('jet')
    # ax.scatter3D(xv, yv, Z, c = Z, cmap = color_map)
    # ax.view_init(90, -90)
    # plt.axis('off')
    # plt.show()
    if type(C1) == str:
        RMS_i = rms_calc(Z)
        C1 = RMS / RMS_i
        Z = C1 * Z
    else:
        Z = C1 * Z
    stat_analysis(Z, N, M, C1, B, sample_type, out_pre)
    return Z


def make_rough(type_sample, z, nodesurf, vertices, angles):
    """
    Applies roughness on the sample.

    :param type_sample: The type of the sample
    :type type_sample: str
    :param z: Surface roughness to apply on the sample
    :type z: array
    :param nodesurf: List of nodes to be moved
    :type nodesurf: array
    :param vertices: Nodes of the sample
    :type vertices: array
    :param angles: List of angles to be followed by roughness (Only in the case of a faceted wire)
    :type angles: array

    :return: Rough sample
    """
    min_dz = abs(z.min())
    if type_sample == 'box' or type_sample == 'grain':
        for i in range(len(z)):
            dz = z[i] + min_dz
            node = nodesurf[i]
            index = node[3]
            poss = np.where(vertices[:, 3] == index)
            vertices[poss, 2] = vertices[poss, 2] + dz
    elif type_sample == 'wire':
        for i in range(len(z)):
            dz = z[i] + min_dz
            node = nodesurf[i]
            thetaa = node[1]
            index = node[3]
            poss = np.where(vertices[:, 3] == index)
            vertices[poss, 0] = vertices[poss, 0] + dz * np.cos(thetaa)
            vertices[poss, 1] = vertices[poss, 1] + dz * np.sin(thetaa)
    elif type_sample == 'poly':
        k = 0
        for i in range(np.shape(z)[0]):
            for j in range(np.shape(z)[1]):
                dz = z[i, j] + min_dz
                node = nodesurf[k]
                index = int(node[3])
                thetaa = np.arctan2(node[1], node[0])
                if thetaa < 0:
                    thetaa = thetaa + 2 * np.pi
                theta_min = np.abs(np.array(angles) - thetaa)
                possi = np.where(abs(theta_min - np.amin(theta_min)) <= 0.01)[0]
                if len(possi) > 1:
                    angle = 0.5*(angles[int(possi[0])]+angles[int(possi[1])])
                    dz = 0.5 * dz
                elif thetaa == 0:
                    angle = 0
                    dz = 0.5 * dz
                else:
                    angle = angles[int(possi)]
                poss = int(np.where(vertices[:, 3] == index)[0])
                vertices[poss, 0] = vertices[poss, 0] + dz * np.cos(angle)
                vertices[poss, 1] = vertices[poss, 1] + dz * np.sin(angle)
                k += 1
    return vertices


def base(radius, nfaces):
    """
    Creates the base of the faceted wire from radius and number of faces.

    :param radius: Radius of the wire
    :type radius: float
    :param nfaces: Number of faces
    :type nfaces: int

    :return: Points forming the base of the wire and angles associated to each face
    """
    points = []
    theta_list = np.linspace(0, 2 * np.pi, nfaces + 1)
    for i in range(len(theta_list)) :
        if theta_list[i] < 0 :
            theta_list[i] += 2*np.pi
    angles = []
    for theta in theta_list[:-1]:
        points.append([radius * np.cos(theta), radius * np.sin(theta)])
    for i in range(len(points)):
        k = i + 1
        if k == len(points):
            k = 0
        x1 = points[int(i)][0]
        x2 = points[int(k)][0]
        y1 = points[int(i)][1]
        y2 = points[int(k)][1]
        the = np.arctan2(0.5 * (y1 + y2), 0.5 * (x1 + x2))
        if the < 0:
            the = the + 2 * np.pi
        angles.append(the)
    angles = np.asarray(angles)
    return points, angles


def cart2cyl(matrix):
    """
    Calculates cylindrical coordinates from cartesian ones.

    :param matrix: List of cartesian coordinates
    :type matrix: array

    :return: Cylindrical coordinates
    """
    cyl_matrix = []
    for p in matrix:
        r = rho(p[0], p[1])
        thet = theta(p[0], p[1])
        cyl_matrix.append([r, thet, p[2], p[3]])
    cyl_matrix = np.asarray(cyl_matrix)
    return cyl_matrix


def cyl2cart(matrix):
    """
    Calculates cartesian coordinates from cylindrical ones.

    :param matrix: List of cylindrical coordinates
    :type matrix: array

    :return: Cartesian coordinates
    """
    cart_matrix = []
    for p in matrix:
        x = p[0] * np.cos(p[1])
        y = p[0] * np.sin(p[1])
        cart_matrix.append([x, y, p[2], p[3]])
    cart_matrix = np.asarray(cart_matrix)
    return cart_matrix


def faces_normals(obj_points, obj_faces):
    """
    Calculates each face's normals.

    :param obj_points: List of points
    :type obj_points: array
    :param obj_faces: List of faces
    :type obj_faces: array

    :return: List of normals
    """
    list_n = []
    for F in obj_faces:
        A = obj_points[int(F[0] - 1)]
        B = obj_points[int(F[1] - 1)]
        C = obj_points[int(F[2] - 1)]
        AB = B - A
        AC = C - A
        n = np.cross(AB, AC)
        n = n / np.linalg.norm(n)
        list_n.append([n[0], n[1], n[2]])
    return list_n


def node_corner(nodesurf):
    """
    From surface nodes, finds all nodes located on edges and corners.

    :param nodesurf: List of surface nodes
    :type nodesurf: array

    :return: List of nodes located on edges and list of nodes located on corners
    """
    all_points = []
    for f in nodesurf:
        for i in f:
            all_points.append(i)
    matrix = np.asarray(all_points)
    list_index = matrix[:, 3]
    value, counts = np.unique(list_index, return_counts=True)
    node_edge = np.where(counts == 2)[0]
    node_corner = np.where(counts >= 3)[0]
    return (node_edge, node_corner)


def make_rough_wulff(vertices, B, C1, RMS, N, M, nodesurf, node_edge, node_corner, list_n):
    """
    Applies roughness on the sample in the case of a Wulff Shaped NP.

    :param vertices: Nodes of the sample
    :type vertices: array
    :param B: The degree of the roughness
    :type B: float
    :param C1: Roughness normalization factor
    :type C1: float
    :param RMS: Root Mean Square
    :type RMS: float
    :param N: Scaling cartesian position
    :type N: int
    :param M: Scaling cartesian position
    :type M: int
    :param nodesurf: List of surface nodes
    :type nodesurf: array
    :param node_edge: List of nodes located on edges
    :type node_edge: array
    :param node_corner: List of nodes located on corners
    :type node_corner: array
    :param list_n: List of face's normals
    :type list_n: list

    :return: Rough Wulff sample
    """
    sfrN, sfrM = vectors(N, M)
    for k in range(len(nodesurf)):
        surf = np.array(nodesurf[k])
        n1 = np.array(list_n[k])
        surf_rot = rot(surf, n1)
        surf_norm = normalize(surf_rot)
        xv = surf_norm[:, 0]
        yv = surf_norm[:, 1]
        m, n = random_numbers(sfrN, sfrM)
        if type(C1) == str:
            z = rdnsurf(m, n, B, xv, yv, sfrM, sfrN)
            C1 = RMS / rms_calc(z)
            z = C1 * z
        else:
            z = C1 * rdnsurf(m, n, B, xv, yv, sfrM, sfrN)
        for i in range(len(surf)):
            p = surf[i]
            index = p[3]
            if index in node_edge:
                delta_z = z[i] / 2
            elif index in node_corner:
                delta_z = z[i] / 3
            else:
                delta_z = z[i]
            poss = np.where(vertices[:, 3] == index)
            vertices[poss, 0] = vertices[poss, 0] + delta_z * n1[0]
            vertices[poss, 1] = vertices[poss, 1] + delta_z * n1[1]
            vertices[poss, 2] = vertices[poss, 2] + delta_z * n1[2]
    return (vertices)


def rot(surf, n1):
    """
    Rotates a surface oriented along n1 axis in order to be oriented along z-axis.

    :param surf: List of nodes from the surface
    :type surf: array
    :param n1: Surface normal
    :type n1: array

    :return: Rotated surface
    """
    n2 = np.array([0, 0, 1])
    n = np.cross(n1, n2)
    if n[0] == 0 and n[1] == 0 and n[2] == 0:
        n = n2
    n = n / np.linalg.norm(n)
    theta = np.arccos(np.dot(n1, n2))
    R = rot_matrix(n, theta)
    surf_rot = []
    for p in surf:
        point = p[:3]
        point_rot = np.dot(R, point)
        surf_rot.append([point_rot[0], point_rot[1], point_rot[2], p[3]])
    surf_rot = np.asarray(surf_rot)
    return (surf_rot)


def rot_matrix(n, theta):
    """
    Generates the rotation matrix. Initial orientation is n, and the angle of rotation is theta. Final orientation is the z-axis.

    :param n: Initial orientation
    :type n: list
    :param theta: Rotation angle
    :type theta: float

    :return: Rotation matrix
    """
    c = np.cos(theta)
    s = np.sin(theta)
    R11 = n[0] * n[0] * (1 - c) + c
    R12 = n[0] * n[1] * (1 - c) - n[2] * s
    R13 = n[0] * n[2] * (1 - c) + n[1] * s
    R21 = n[0] * n[1] * (1 - c) + n[2] * s
    R22 = n[1] * n[1] * (1 - c) + c
    R23 = n[1] * n[2] * (1 - c) - n[0] * s
    R31 = n[0] * n[2] * (1 - c) - n[1] * s
    R32 = n[1] * n[2] * (1 - c) + n[0] * s
    R33 = n[2] * n[2] * (1 - c) + c
    R = np.array([[R11, R12, R13], [R21, R22, R23], [R31, R32, R33]])
    return (R)


def normalize(surf):
    """
    Normalizes the coordinates of points composing the surface.

    :param surf: List of nodes of the surface
    :type surf: array

    :return: Normalized surface
    """
    X = surf[:, 0]
    Y = surf[:, 1]
    Z = surf[:, 2]
    T = surf[:, 3]
    x_max = np.max(abs(X))
    y_max = np.max(abs(Y))
    Xf = (X / x_max + 1) / 2
    Yf = (Y / y_max + 1) / 2
    Zf = 0 * Z
    surf_norm = [[Xf[i], Yf[i], Zf[i], T[i]] for i in range(len(Xf))]
    surf_norm = np.asarray(surf_norm)
    return (surf_norm)


def cube_faces(length):
    """
    Generates the points and faces of a cube

    :param length: Size of the cube
    :type length: float

    :return: Points and faces of the cube
    """
    obj_points = [[0, 0, 0],
                  [length, 0, 0],
                  [0, length, 0],
                  [length, length, 0],
                  [0, 0, length],
                  [length, 0, length],
                  [0, length, length],
                  [length, length, length]]
    obj_faces = [[1, 2, 3, 4],
                 [1, 2, 5, 6],
                 [2, 4, 6, 8],
                 [4, 8, 3, 7],
                 [1, 5, 3, 7],
                 [5, 6, 7, 8]]
    obj_faces = np.asarray(obj_faces)
    obj_points = np.asarray(obj_points)
    return obj_points, obj_faces


def rotate_obj_wulff(obj_points, orien_x, orien_z):
    """
    :param obj_points: Points describing the facets of the wulff-shape
    :type obj_points: array
    :param orien_x: Orientation along x-axis
    :type orien_x: list
    :param orien_z: Orientation along z-axis
    :type orien_z: list

    :return: Points respecting the desired orientation
    """
    n2 = np.array([0, 0, 1])
    n = np.cross(orien_z, n2)
    if n[0] == 0 and n[1] == 0 and n[2] == 0:
        n = n2
    n = n / np.linalg.norm(n)
    theta = np.arccos(np.dot(orien_z, n2) / (np.linalg.norm(orien_z) * np.linalg.norm(n2)))
    R = rot_matrix(n, theta)
    surf_rot = []
    for p in obj_points:
        point_rot = np.dot(R, p)
        surf_rot.append([point_rot[0], point_rot[1], point_rot[2]])
    surf_rot = np.asarray(surf_rot)
    x_axis = np.array([1, 0, 0])
    x_rot = np.dot(R, x_axis)
    theta_x = np.arccos(np.dot(x_rot, orien_x) / (np.linalg.norm(x_rot) * np.linalg.norm(orien_x)))
    R_x = np.array([[np.cos(theta_x), -1 * np.sin(theta_x), 0], [np.sin(theta_x), np.cos(theta_x), 0], [0, 0, 1]])
    points_f = []
    for p in surf_rot:
        rot_f = np.dot(R_x, p)
        points_f.append([rot_f[0], rot_f[1], rot_f[2]])
    obj_points_f = np.asarray(points_f)
    # fig = plt.figure(figsize=(4, 4))
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(obj_points_f[:,0], obj_points_f[:,1], obj_points_f[:,2])
    # plt.show()
    return (obj_points_f)

def align_poly(vertices, angles):
    """
    Aligns the first facet with x-axis

    :param vertices: List of nodes
    :type vertices: array
    :param angles: List of angles
    :type angles: array

    :return: Rotated nodes
    """
    theta = angles[0]
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    rotated_points = np.dot(vertices, R)
    return (rotated_points)

####################################################################
####################### SURFACE ANALYSIS ###########################
####################################################################

def rough(x, y, H, C1, M, N):
    """
    Generates a rough surface

    :param x: List of x coordinates
    :type x: array
    :param y: List of y coordinates
    :type y: array
    :param H: Hurst exponent
    :type H: float
    :param C1: Normalization factor
    :type C1: float
    :param M: Scaling cartesian position
    :type M: int
    :param N: Scaling cartesian position
    :type N: int

    :return: Height matrix
    """
    z = np.zeros(np.shape(x))
    listM = np.linspace(-M, M, 2 * M + 1)
    listN = np.linspace(-N, N, 2 * N + 1)
    k = 0
    for m in listM:
        for n in listN:
            if m == 0 and n == 0:
                pass
            else:
                mod = (m * m + n * n) ** (-1 * (1 + H))
                G = np.random.randn()
                U = -np.pi / 2 + np.pi * np.random.rand()
                zadd = C1 * G * mod * np.cos(2 * np.pi * (m * x + n * y) + U)
                z = z + zadd
            k += 1
    return (z)


def sigma(z):
    """
    Calculates the standard deviation of a height distribution

    :param z: height matrix
    :type z: array

    :return: Standard deviation
    """
    Z = z.flatten()
    l = np.mean(Z)
    return (round(np.sqrt(np.sum(np.power(Z - l, 2)) / len(Z)), 3))


def sk(z):
    """
    Calculates the skewness of a height distribution

    :param z: height matrix
    :type z: array

    :return: Skewness
    """
    Z = z.flatten()
    l = np.mean(Z)
    return (round((1 / (sigma(z) ** 3 * len(Z))) * np.sum(np.power(Z - l, 3)), 3))


def Kurto(z):
    """
    Calculates the kurtosis of a height distribution

    :param z: height matrix
    :type z: array

    :return: Kurtosis
    """
    Z = z.flatten()
    l = np.mean(Z)
    return (round((1 / (sigma(z) ** 4 * len(Z))) * np.sum(np.power(Z - l, 4)), 3))


def rescale(D, scale):
    """
    Rescales the height distribution between zmin and zmax

    :param D: Height matrix
    :type D: array
    :param scale: zmin and zmax
    :type scale: list

    :return: Rescaled matrix
    """
    D = (D - np.min(D))
    D = D / np.max(D)
    lower = scale[0]
    upper = scale[1]
    Df = [lower + (upper - lower) * x for x in D]
    return (np.asarray(Df))

def box_perio(file_lmp, dis):
    fint = open(file_lmp, "r")
    lines = fint.readlines()
    data_x = lines[5].split( )
    lines[5] = '{} {} {} {}\n'.format(data_x[0], str(float(data_x[1])+dis[0]), data_x[2], data_x[3])
    data_y = lines[6].split( )
    lines[6] = '{} {} {} {}\n'.format(data_y[0], str(float(data_y[1])+dis[1]), data_y[2], data_y[3])
    data_z = lines[7].split( )
    lines[7] = '{} {} {} {}\n'.format(data_z[0], str(float(data_z[1])+dis[2]), data_z[2], data_z[3])
    fint.close()
    fend = open(file_lmp, "w")
    fend.writelines(lines)
    fend.close()

def refine(stl, ext):
    print('====== > Refining the Mesh')
    filee = Path(stl).stem
    mesh = meshio.read(stl)
    vertices, faces = mesh.points, mesh.cells
    faces = faces[0][1]
    with pygmsh.geo.Geometry() as geom:
        for f in faces:
            poly = geom.add_polygon([vertices[f[0]], vertices[f[1]], vertices[f[2]]], mesh_size=3)
        mesh = geom.generate_mesh()
        vertices = mesh.points
        faces = mesh.get_cells_type('triangle')
    # ext = ['obj', 'stl', 'inp', 'vtk']
    for e in ext:
        meshio.write_points_cells(filee+'.'+e, vertices, {"triangle": faces})