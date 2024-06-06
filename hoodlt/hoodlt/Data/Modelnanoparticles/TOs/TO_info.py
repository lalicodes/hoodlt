"""
:module: TO_info
:platform: Unix, Windows
:synopsis: Calculates the radius of the given truncated octahedron (TO),
gives a list of available TO's within a range of radius

.. moduleauthor:: Xun Zha <xzha@iastate.edu> November 2018
"""

import numpy as np
import numpy.linalg as la
import pkg_resources
import pandas as pd
import hoodlt.Utils.neighbor_number as nn


class ConfigInfo(object):
    def __init__(self, num_pnts, version, atom_type='Au'):
        """

        :param num_pnts:
        :param version:
        :param atom_type:
        """

        # lattice constant
        filename = pkg_resources.resource_filename('hoodlt', 'Data/Modelnanoparticles/TOs/fcc_lattice_constant.csv')
        l_data = pd.read_csv(filename)
        l_constant = l_data[l_data.loc[:, 'atom'] == atom_type]['lattice constant'].values[0]
        self.a_nn = l_constant/np.sqrt(2)  # nearest neighbor distance

        # filename with the TO information
        name = 'N'+str(int(num_pnts))+'_v'+str(int(version))+'.txt'
        filename = pkg_resources.resource_filename('hoodlt', 'Data/Modelnanoparticles/TOs/'+name)

        # read atom coordinates
        pos = np.genfromtxt(filename)  # all atom coordinates
        n_list = nn.neighbor_number(pos, 1)  # list of coordination number for each atom
        position = pos[n_list < 12]  # coordinates of surface atoms

        # vectors normal to (100) and (111) facets
        vec100 = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        vec111 = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1],
                           [1, 1, -1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])/np.sqrt(3)
        rotation45 = np.array([[1/np.sqrt(2), -1/np.sqrt(2), 0],
                               [1/np.sqrt(2), 1/np.sqrt(2), 0],
                               [0, 0, 1]])  # NCs were rotated by 45 degrees around z-axis
        vec100 = np.matmul(rotation45, vec100.T).T
        vec111 = np.matmul(rotation45, vec111.T).T

        # list of atoms coordinates and numbers on the six (100) facets
        temp = np.matmul(vec100, position.T)
        dis100 = np.amax(temp, axis=1)
        arg100 = np.array([np.isclose(temp[i], dis100[i]) for i in range(len(dis100))])
        pos100 = list([position[arg100[i]] for i in range(len(dis100))])
        # area of (100) facets
        num100 = np.sum(arg100, axis=1)  # mxn
        l100 = np.zeros(len(dis100))  # diagonal length of (100) rectangles
        for i, pos in enumerate(pos100):
            mat = np.zeros([len(pos), len(pos), 3])
            for j in range(len(pos)):
                mat[j] = pos - pos[j]
            dis = la.norm(mat, axis=-1)
            l100[i] = np.amax(dis)  # sqrt((m-1)^2+(n-1)^2)
        area100 = np.sum(num100-np.sqrt(l100**2+2*num100-1))  # (m-1)x(n-1) = mxn-(m+n-1) = mxn-sqrt(l^2+2mxn-1)
        # list of atom coordinates and numbers on the eight (111) facets
        temp = np.matmul(vec111, position.T)
        dis111 = np.amax(temp, axis=1)
        arg111 = np.array([np.isclose(temp[i], dis111[i]) for i in range(len(dis111))])
        pos111 = list([position[arg111[i]] for i in range(len(dis111))])
        # area of (111) facets
        area111 = 0
        for i, pos in enumerate(pos111):
            temp = nn.neighbor_number(pos, 1)
            corner = sum(temp == 3)
            side = sum(temp == 4)
            inner = len(pos) - corner - side
            area111 += (corner/3. + side/2. + inner)*np.sqrt(3)/2.
        # surface area
        self.area = area100 + area111
        # radius from surface area
        self.radius = np.sqrt(self.area/4./np.pi)
        # volume
        corner = sum(n_list == 6)
        side = sum(n_list == 7)
        surface = sum(n_list == 8) + sum(n_list == 9)
        inner = sum(n_list == 12)
        self.volume = self.a_nn**3/np.sqrt(2) * (corner/4.+side/3.+surface/2.+inner)
        # radius from volume
        self.radius_v = np.cbrt(self.volume*3/4./np.pi)


def find_n(radius, atom_type='Au'):
    """find TOs with given radius or given radius range

    :param radius:
    :param atom_type:
    :return: a dictionary of TOs' number, version, radius information
    """
    # lattice constant
    filename = pkg_resources.resource_filename('hoodlt', 'Data/Modelnanoparticles/TOs/fcc_lattice_constant.csv')
    l_data = pd.read_csv(filename)
    l_constant = l_data[l_data.loc[:, 'atom'] == atom_type]['lattice constant'].values[0]

    # read list of available TOs
    filename = pkg_resources.resource_filename('hoodlt', 'Data/Modelnanoparticles/TOs/number_radius.txt')
    data = np.genfromtxt(filename)
    data[:, 2] *= l_constant/np.sqrt(2)

    # sort list of TOs by radius, number, version
    arg = np.lexsort((data[:, 1], data[:, 0], data[:, 2]))
    data = data[arg]

    # select the TO with given radius or TOs within the radius range
    if isinstance(radius, float) or isinstance(radius, int):
        arg = np.abs(data[:, 2]-radius).argmin()
        number, version, radius = data[arg]
        dict_info = {'number': int(number), 'version': int(version), 'radius': radius}
    elif len(radius) == 2:
        arg = np.all([data[:, 2] > radius[0], data[:, 2] < radius[1]], axis=0)
        number, version, radius = data[arg].T
        dict_info = {}
        for i in range(len(number)):
            dict_info[str(i)] = {'number': int(number[i]), 'version': int(version[i]), 'radius': radius[i]}
    else:
        raise ValueError('Input radius needs to be a number, or a range, check input radius')

    return dict_info
