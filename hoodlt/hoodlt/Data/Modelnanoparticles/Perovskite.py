"""
:module: Perovskite
:platform: Unix, Windows
:synopsis: Defines the class for perovskite structure, ABX3

.. moduleauthor:: Xun Zha <xzha@iastate.edu> January 2018
.. history:
..              Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
..              Xun Zha <xzha@iastate.edu> Sept. 2020
..                  - minor corrections
..              Xun Zha <xzha@iastate.edu> Dec. 2020
..                  - corrections
..              Xun Zha <xzha@iastate.edu> June 2021
..                  - removed parameter ligand_name, added parameter ligand_position
..              Xun Zha <xzha@iastate.edu> June 2021
..                  - Reading of the file consistent with the units in the forcefield (adding Angstrom to Construction)
"""

import numpy as np
import copy
from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs


class Perovskite(NanoAbs):
    """
    Defines the Perovskite Nanocrystal

    if :math:`l` is the linear size there are :math:`l^3` Wigner Seitz unit cells.

    There are :math:'(l+1)^3` Cs, math:`l^3` Pb and math:`3(l+1)l^2` Br atoms in total (or equivalent)

    Therefore there are a maximum of :math:`6l^2+2` Cs and :math:`6l^2` Br binding sides (or equivalent)
    """
    def __init__(self, forcefield, l_size, atom_type=None, graft_type=None, ligand_position=None,
                 core_name_suffix=None):
        """

        :param forcefield: name of the forcefield
        :param l_size: linear size (linear number of Perovskite unit cells)
        :param atom_type: atom types of the NC core (default are Cs, Pb, Br)
        :param graft_type: linkers as (A,B) A links to Cs, B links to Br
        :param ligand_position: list of two numpy arrays of the indices for the positions that bind to Cs, and
                                the indices for the positions that bind to Br
        :param core_name_suffix: a suffix added to the end of the core name when there are different rigid cores
        """

        # core atom types
        if atom_type is None:
            atom_type = ['Cs', 'Pb', 'Br']  # default core atom types
        elif len(atom_type) != 3:
            raise ValueError('Check atom type input')

        # grafting atom type
        if graft_type is None:
            graft_type = ['Coo', 'NH3']  # default graft atom types
        elif len(graft_type) != 2:
            raise ValueError('Check graft_type input, needs to be a list of two items')

        # maximum number of grafted atoms
        graft_max_num = np.array([6 * l_size ** 2 + 2, 6 * l_size ** 2])
        # actual number of grafted atoms
        if ligand_position is None:
            graft_num = graft_max_num
        elif len(ligand_position) == 2:
            if len(ligand_position[0]) > graft_max_num[0] or len(ligand_position[1]) > graft_max_num[1]:
                raise ValueError('Number of indices for ligand positions exceed the maximum possible grafting sites')
            graft_num = np.array([len(ligand_position[0]), len(ligand_position[1])])
        else:
            raise ValueError('Parameter ligand_position needs to be a list of two numpy arrays')

        # core particles
        self._core_total_num = np.array([np.power(l_size+1, 3), np.power(l_size, 3), 3*(l_size+1)*(l_size**2)])
        self._core_surface_num = np.array([6 * l_size ** 2 + 2, 0, 6 * l_size ** 2])

        # number of core atoms (rigid center plus constituent particles)
        num_particles = 1 + np.sum(self._core_surface_num)
        # nano name
        core_name = ''.join(atom_type) + '3' + graft_type[0] + str(graft_num[0]) + graft_type[1] + str(graft_num[1])
        if core_name_suffix is not None:
            if not isinstance(core_name_suffix, (str, int)):
                raise ValueError('Parameter "core_name_suffix" needs to be a string or an integer')
            core_name += '_' + str(core_name_suffix)

        super(Perovskite, self).__init__(forcefield, num_particles, core_name)

        # core types
        self.types = ['_'+self.name] + [atom_type[0]] + [atom_type[2]]
        self.typeid = ['_'+self.name]+[atom_type[0]]*self._core_surface_num[0]+[atom_type[2]]*self._core_surface_num[2]

        # graft info
        self.graft_type = graft_type
        self.graft_num = graft_num

        # lattice size
        self._l_size = l_size
        # lattice constant
        self.params_outside_forcefield['Cs-Cs'] = 5.87
        self._l_const = self.params_outside_forcefield['Cs-Cs'] = 5.87
        # get the units of the forcefield
        angstrom_to_construction = self.ff_reader.get_units().angstrom_to_construction

        # positions of the core atoms
        pos_all, pos_surf = self._core_position(self._l_size, self._l_const)
        self.position = np.vstack((np.zeros([1, 3]), pos_surf)) * angstrom_to_construction

        # positions of the grafted atoms
        self.graft_sites = self._grafting_position(copy.copy(pos_surf), ligand_position) * angstrom_to_construction

        # mass and moment of inertia
        mass_all, mass_surf = self._mass(atom_type)
        self.mass = np.append([np.sum(mass_all)], mass_surf)
        self.moment_inertia[0] = np.diag(self.moment_of_inertia())

    @staticmethod
    def _core_position(l_size, l_const):
        """calculate the positions of core atoms

        :param l_size: lattice size
        :param l_const: lattice constant
        :return: lists of positions of core atoms
        """

        # basis vectors
        v_vector = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])

        # positions of all atoms in the NC core
        pos_a = []
        pos_b = []
        pos_x = []
        for ind_1 in range(l_size):
            for ind_2 in range(l_size):
                for ind_3 in range(l_size):
                    ary = np.array([ind_1, ind_2, ind_3])
                    pos_a.append(v_vector[0]+ary)
                    pos_b.append(v_vector[1]+ary)
                    pos_x.append(v_vector[2:]+ary)
                pos_a.append(np.array([[ind_1, ind_2, l_size], [ind_1, l_size, ind_2], [l_size, ind_1, ind_2]]))
                pos_x.append(np.array([[ind_1+0.5, ind_2+0.5, l_size], [ind_1+0.5, l_size, ind_2+0.5],
                                       [l_size, ind_1+0.5, ind_2+0.5]]))
            pos_a.append(np.array([[ind_1, l_size, l_size], [l_size, ind_1, l_size], [l_size, l_size, ind_1]]))
        pos_a.append(l_size*np.ones(3))

        # positions of surface atoms in the NC core (constituent particles of the rigid core)
        pos_a_s = []
        pos_x_s = []
        for ind_1 in range(l_size):
            for ind_2 in range(l_size):
                pos_a_s.append(np.array([[ind_1, ind_2, 0], [ind_1+1, ind_2+1, l_size],
                                         [ind_1+1, 0, ind_2+1], [ind_1, l_size, ind_2],
                                         [0, ind_1, ind_2+1], [l_size, ind_1+1, ind_2]]))
                pos_x_s.append(np.array([[ind_1+0.5, ind_2+0.5, 0], [ind_1+0.5, ind_2+0.5, l_size],
                                         [ind_1+0.5, 0, ind_2+0.5], [ind_1+0.5, l_size, ind_2+0.5],
                                         [0, ind_1+0.5, ind_2+0.5], [l_size, ind_1+0.5, ind_2+0.5]]))
        pos_a_s.append(np.array([[l_size, 0, 0], [0, l_size, l_size]]))

        pos_all = np.vstack((np.vstack(pos_a), pos_b, np.vstack(pos_x)))*l_const
        pos_all -= np.mean(pos_all, axis=0)
        pos_surf = np.vstack((np.vstack(pos_a_s), np.vstack(pos_x_s)))*l_const
        pos_surf -= np.mean(pos_surf, axis=0)

        return pos_all, pos_surf

    def _grafting_position(self, pos, ligand_position):
        """calculate the position of grafting atoms

        :param pos: positions of surface core atoms
        :param ligand_position: a list of two arrays of indices for the grafted positions
        :return: positions of the graft atoms
        """
        # Cs(A)-Coo(graft1), Br(X)-NH3(graft2) bond length
        self.params_outside_forcefield['Cs-Coo'] = 5.08356912021466
        self.params_outside_forcefield['Br-NH3'] = 2.55001491346532
        angstrom_to_construction = self.ff_reader.get_units().angstrom_to_construction
        length_a = self.params_outside_forcefield['Cs-Coo']*angstrom_to_construction
        length_x = self.params_outside_forcefield['Br-NH3']*angstrom_to_construction

        # Cs(A) and Br(X) positions
        pos_a = pos[:self._core_surface_num[0]]
        pos_x = pos[self._core_surface_num[0]:]
        if ligand_position is not None:
            pos_a = pos_a[ligand_position[0]]
            pos_x = pos_x[ligand_position[1]]

        # normal vectors of the six facets
        vec = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])

        # distance of the six (100) facets to the origin
        dis = np.amax(np.matmul(vec, pos.T), axis=1)

        # find Cs's and calculate Coo's positions
        temp = np.matmul(vec, pos_a.T)  # Cs's at the six (100) facets
        arg = np.array([np.isclose(temp[i], dis[i]) for i in range(len(dis))])  # lists of Cs's at the six (100) facets
        corner_count = np.sum(arg, axis=0)  # Cs at the corners count 3, at the sides count 2, at the surface count 1
        for i in range(len(vec)):
            pos_a[arg[i]] += np.outer(1/np.sqrt(corner_count[arg[i]]), vec[i]) * length_a

        # find Br's and calculate NH3's positions
        temp = np.matmul(vec, pos_x.T)  # Br's at the six (100) facets
        arg = np.array([np.isclose(temp[i], dis[i]) for i in range(len(dis))])  # lists of Br's at the six (100) facets
        for i in range(len(vec)):
            pos_x[arg[i]] += vec[i] * length_x

        return np.vstack((pos_a, pos_x))

    def _mass(self, atom_type):
        """

        :return: lists of mass of core atoms
        """
        mass_a = self.ff_reader.get_molecular_weight(atom_type[0])
        mass_b = self.ff_reader.get_molecular_weight(atom_type[1])
        mass_x = self.ff_reader.get_molecular_weight(atom_type[2])

        # calculate the mass of the core bulk atoms
        mass_all = np.concatenate((np.ones((self._core_total_num[0], 1))*mass_a,
                                   np.ones((self._core_total_num[1], 1))*mass_b,
                                   np.ones((self._core_total_num[2], 1))*mass_x))

        # calculate the mass of the core surface atoms
        mass_surface = np.concatenate((np.ones((self._core_surface_num[0], 1))*mass_a,
                                       np.ones((self._core_surface_num[2], 1))*mass_x))

        return mass_all, mass_surface

    def volume(self):
        """

        :return: value of the volume
        """
        return np.power(self._l_size*self._l_const, 3)

    def area(self):
        """

        :return: value of the area
        """
        return np.power(self._l_size*self._l_const, 2)
