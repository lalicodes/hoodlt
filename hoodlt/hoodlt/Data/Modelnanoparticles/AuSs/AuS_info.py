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
    def __init__(self, num_pnts, version, atom_type):
        """

        :param num_pnts: number of atoms in the TO configuration
        :param version: degeneracy
        :param atom_type: type of atom in the NC
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
        pos_surf = pos[n_list < 12]  # coordinates of surface atoms

        # volume
        s_edge = sum(n_list == 4)  # single atom on a (100) facet of a small NC
        e_edge = sum(n_list == 5)  # edge atom on a (100) or (111) facet of a small NC, double check ???
        edge = sum(n_list == 6)  # atoms at the edge with 6 neighbors
        side = sum(n_list == 7)  # atoms at the side with 7 neighbors
        surface = sum(n_list == 8) + sum(n_list == 9)  # atoms at the surface with 8 or 9 neighbors
        inner = sum(n_list == 12)  # atoms in the bulk with 12 neighbors
        self.volume = (s_edge/9.+e_edge*7/36.+edge/4.+side/3.+surface/2.+inner)/np.sqrt(2)*self.a_nn**3
        # radius from volume
        self.radius_v = np.cbrt(self.volume*3/4./np.pi)

        # vectors normal to (100) and (111) facets
        self.vec100 = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        self.vec111 = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1],
                                [1, 1, -1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])/np.sqrt(3)
        # rotate the vectors by 45 degrees around z-axis
        rotation45 = np.array([[1/np.sqrt(2), -1/np.sqrt(2), 0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0, 0, 1]])
        self.vec100 = np.matmul(rotation45, self.vec100.T).T
        self.vec111 = np.matmul(rotation45, self.vec111.T).T

        # list of atom coordinates on the six (100) facets, area of the six (100) facets
        temp = np.matmul(self.vec100, pos_surf.T)
        dis100 = np.amax(temp, axis=1)
        arg100 = np.array([np.isclose(temp[i], dis100[i]) for i in range(len(dis100))])
        pos100 = list([pos_surf[arg100[i]] for i in range(len(dis100))])
        # area of (100) facets
        n_temp = list([nn.neighbor_number(pos100[i], 1) for i in range(len(dis100))])
        corner100 = np.array([sum(n_temp[i] == 2) for i in range(len(dis100))])
        side100 = np.array([sum(n_temp[i] == 3) for i in range(len(dis100))])
        inner100 = np.array([sum(n_temp[i] == 4) for i in range(len(dis100))])
        n_squares = corner100/4.+side100/2.+inner100
        # list of atom coordinates on the eight (111) facets, area of the eight (111) facets
        temp = np.matmul(self.vec111, pos_surf.T)
        dis111 = np.amax(temp, axis=1)
        arg111 = np.array([np.isclose(temp[i], dis111[i]) for i in range(len(dis111))])
        pos111 = list([pos_surf[arg111[i]] for i in range(len(dis111))])
        # area of (111) facets
        n_temp = list([nn.neighbor_number(pos111[i], 1) for i in range(len(dis111))])
        t_corner111 = np.array([sum(n_temp[i] == 2) for i in range(len(dis111))])  # corner of a triangle
        corner111 = np.array([sum(n_temp[i] == 3) for i in range(len(dis111))])
        side111 = np.array([sum(n_temp[i] == 4) for i in range(len(dis111))])
        inner111 = np.array([sum(n_temp[i] == 6) for i in range(len(dis111))])
        n_triangles = 2*(t_corner111/6.+corner111/3.+side111/2+inner111)

        # surface area
        self.area = (np.sum(n_squares) + np.sum(n_triangles)*np.sqrt(3)/4.)*self.a_nn**2
        # radius from surface area
        self.radius = np.sqrt(self.area/4./np.pi)
        # calculate positions of grafting atoms
        # get forcefield parameters
        force_field = pkg_resources.resource_filename('hoodlt', 'Data/Forcefield/params_forcefield.xlsx')
        xls_ff = pd.ExcelFile(force_field)
        sig_s = xls_ff.parse('nonbonded')[xls_ff.parse('nonbonded').loc[:, 'name'] == atom_type+'-S']['sigma'].values[0]
        sig_s2 = xls_ff.parse('nonbonded')[xls_ff.parse('nonbonded').loc[:, 'name'] == 'S-S']['sigma'].values[0]
        dis_s111 = np.sqrt(sig_s**2-self.a_nn**2/3)  # ???
        dis_s100 = np.sqrt(sig_s**2-self.a_nn**2/2)  # ???
        dis_s_s = np.sqrt(sig_s**2-self.a_nn**2/4)  # ???
        dis_s_e = sig_s  # ???
        # print(dis_s100, dis_s111, dis_s_s, dis_s_e)
        area100_p = n_squares * (dis100+dis_s100/self.a_nn)/dis100 * self.a_nn**2
        area111_p = n_triangles*np.sqrt(3)/4. * (dis111+dis_s111/self.a_nn)/dis111 * self.a_nn**2
        # surface area of sulfurs
        self.area_sulfur = np.sum(area100_p)+np.sum(area111_p)

        '''
        # list of atom coordinates on the edges and sides
        pos_edge = pos[n_list == 6]
        pos_side = pos[np.any([n_list == 6, n_list == 7], axis=0)]

        # hollow sites on edges
        dis_eg = la.norm(pos_edge, axis=-1)*self.a_nn
        pos_graft = np.array([pos_edge[i]*(dis_eg[i]+dis_s_e)/dis_eg[i] for i in range(len(dis_eg))])*self.a_nn
        index_g = np.repeat(0, len(pos_edge))

        # hollow sites on sides
        graft_side = []
        for ind_1 in range(len(pos_side)):
            dis = la.norm(pos_side-pos_side[ind_1], axis=-1)
            arg = np.where(np.isclose(dis, 1))[0]
            for ind_2 in arg:
                if ind_2 > ind_1:
                    graft_side.append(np.mean([pos_side[ind_1], pos_side[ind_2]], axis=0))
        dis_sd = la.norm(graft_side, axis=-1)*self.a_nn
        pos_graft = np.vstack((pos_graft, np.array([graft_side[i]*(dis_sd[i]+dis_s_s)/dis_sd[i]
                                                    for i in range(len(dis_sd))])*self.a_nn))
        index_g = np.append(index_g, np.repeat(1, len(dis_sd)))

        # hollow site on (100) facets
        graft100 = list([np.zeros([int(n_squares[i]), 3]) for i in range(len(n_squares))])
        for i, pos in enumerate(pos100):
            temp = []
            for ind_1 in range(len(pos)):
                dis = la.norm(pos-pos[ind_1], axis=-1)
                arg10 = np.isclose(dis, 1)
                neighbor = np.where(np.isclose(dis, np.sqrt(2)))[0]
                for ind_2 in neighbor:
                    dis = la.norm(pos-pos[ind_2], axis=-1)
                    arg20 = np.isclose(dis, 1)
                    arg = np.where(np.all([arg10, arg20], axis=0))[0]
                    temp.append([ind_1, ind_2, arg[0], arg[1]])
            temp_list = np.unique([sorted(temp[i]) for i in range(len(temp))], axis=0)
            if len(temp_list) != int(n_squares[i]):
                raise ValueError('Check (100) squares counting!')
            for j in range(len(temp_list)):
                arg = temp_list[j]
                graft100[i][j] = np.mean([pos[arg[0]], pos[arg[1]], pos[arg[2]], pos[arg[3]]], axis=0)
        for i in range(len(graft100)):
            pos_graft = np.vstack((pos_graft, graft100[i]*self.a_nn + dis_s100*self.vec100[i]))
            index_g = np.append(index_g, np.repeat(2, len(graft100[i])))
        for i, pos in enumerate(pos100):
            pos = pos[nn.neighbor_number(pos, 1) > 2]
            pos_graft = np.vstack((pos_graft, pos*self.a_nn + dis_s_e*self.vec100[i]))
            index_g = np.append(index_g, np.repeat(2, len(pos)))

        # hollow sites on (111) facets
        graft111 = list([np.zeros([int(n_triangles[i]), 3]) for i in range(len(n_triangles))])
        for i, pos in enumerate(pos111):
            temp = []
            for ind_1 in range(len(pos)):
                dis = la.norm(pos-pos[ind_1], axis=-1)
                arg = np.isclose(dis, 1)
                neighbor_1 = np.where(arg)[0]
                for ind_2 in neighbor_1:
                    dis = la.norm(pos-pos[ind_2], axis=-1)
                    arg = np.isclose(dis, 1)
                    neighbor_2 = np.where(arg)[0]
                    for ind_3 in neighbor_2:
                        if np.isin(ind_3, neighbor_1):
                            temp.append([ind_1, ind_2, ind_3])
            temp_list = np.unique([sorted(temp[i]) for i in range(len(temp))], axis=0)
            if len(temp_list) != int(n_triangles[i]):
                raise ValueError('Check (111) triangles counting!')
            for j in range(len(temp_list)):
                arg = temp_list[j]
                graft111[i][j] = np.mean([pos[arg[0]], pos[arg[1]], pos[arg[2]]], axis=0)
        for i in range(len(graft111)):
            pos_graft = np.vstack((pos_graft, graft111[i]*self.a_nn + dis_s111*self.vec111[i]))
            index_g = np.append(index_g, np.repeat(3, len(graft111[i])))

        dis = np.zeros([len(pos_graft), len(pos_graft)])
        for i in range(len(pos_graft)):
            dis[i] = la.norm(pos_graft-pos_graft[i], axis=-1)

        # criteria = np.array([1, 1, np.sqrt(5/2), np.sqrt(3)])*sig_s2
        criteria = np.ones(4)*sig_s2
        n1 = 0
        n2 = np.where(dis[n1] >= criteria[index_g[n1]])[0][0]
        l_trig = [n1, n2]
        for i in range(len(pos_graft)):
            arg = np.where(np.all([dis[n1] >= criteria[index_g[n1]], dis[n2] >= criteria[index_g[n2]]], axis=0))[0]
            temp = []
            for x in arg:
                if np.all(dis[l_trig, x] >= criteria[0]):
                    temp.append(x)
            if len(temp) == 0:
                break
            arg = np.lexsort((dis[n2, temp], dis[n1, temp]))
            n1 = n2
            n2 = temp[arg[0]]
            l_trig.append(n2)
        '''

        num_graft = np.ceil(self.area_sulfur/10).astype(int)
        pos_graft = np.zeros((num_graft, 3))

        phi = np.pi*(3-np.sqrt(5))*np.arange(num_graft)
        cs = (2*np.arange(num_graft)+1)/num_graft-1
        sn = np.sqrt(1-cs**2)

        pos_graft[:, 0] = sn*np.cos(phi)
        pos_graft[:, 1] = sn*np.sin(phi)
        pos_graft[:, 2] = cs

        max_dis = np.amax(la.norm(pos_surf, axis=1)*self.a_nn)
        pos_graft *= max_dis+4

        self.pos_all = pos*self.a_nn
        self.pos_surf = pos_surf*self.a_nn
        self.pos_graft = pos_graft

        # self.pos_graft = pos_graft[l_trig]
        self.pos100 = [pos*self.a_nn for pos in pos100]
        self.pos111 = [pos*self.a_nn for pos in pos111]
        self.dis100 = dis100
        self.dis111 = dis111


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
