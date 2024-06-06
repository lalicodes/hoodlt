"""
:module: PEO_bis
:platform: Unix, Windows
:synopsis: Implements a class defining polyethylene oxide

.. moduleauthor:: Elizabeth Macias <emacias@iastate.edu> March 2022
.. history:
..                Elizabeth Macias <emacias@iastate.edu> March 2022
..                  - center of mass is outside build_chain() function
..                  - all CH3 angles are ~109
..                  - more than 2-units of PEO
..
..                Alex Travesset
..                  - simplified the class
..
..                Elizabeth Macias <emacias@iastate.edu> September 2022
..                  - added sulfur end to PEO to make it a ligand
"""


import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class PEOLigand(LigandAbs):
    """
    Defines a PEO polymer
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        # there are 2 carbons 1 oxygen and 4 hydrogens per monomer
        self.repeat_length = 7

        super(PEOLigand, self).__init__(repeats, ff, repeats*self.repeat_length+2+(3), 'PEOLigand')

        '''
        7 particles in one repeat unit, 4 particles in the ch3 rigid body, and 1 oxygen outside of the repeat-chain
        '''

        self.diameter = np.zeros(self.num_particles) + 0.08
        self.types = ['C2', 'C3', 'HP', 'OP', 'CS2', 'HA', 'SD']
        self.end_group = 2 * [None]
        self.end_group[0] = ['SD', 'HA', 'HA']+ ['CS2', 'HP', 'HP', 'C2'] + ['OP', 'C2', 'HP', 'HP']
        self.end_group[1] = ['HP', 'HP', 'C2', 'OP', 'C3', 'HP', 'HP', 'HP']
        if repeats == 1:
            self.typeid = ['SD', 'HA', 'HA', 'CS2', 'HP', 'HP', 'C2', 'OP', 'C3', 'HP', 'HP', 'HP']
        else:
            self.typeid = self.end_group[0] + (repeats-2)*['HP', 'HP', 'C2', 'OP', 'C2', 'HP', 'HP'] + self.end_group[1]

        self.mass = np.zeros(self.num_particles)
        self.charge = np.zeros(self.num_particles)
        self.body = self.num_particles*[-1]

        for ind, typ in enumerate(self.typeid):
            self.mass[ind] = self.ff_reader.get_molecular_weight(typ)
            self.charge[ind] = self.ff_reader.get_charge(typ)

        # particle positions
        self.position = np.zeros([self.num_particles-3, 3])
        self.build_chain()


        new = np.zeros([self.num_particles, 3])
        new[4:] = self.position[1:]
        new[:4] = self.position[:4] - (self.position[5]-self.position[4])
        zrotation = np.array(( (-1, 0, 0),(0, -1, 0), (0, 0, 1) ))
        new[:4] = [zrotation.dot(p) for p in new[:4]]

        self.position = new

    def build_chain(self):
        """
        called by initializer to get the positions
        
        :return: positions relative to first C
        """

        r_c3h = self.ff_reader.get_bond_r0('C3-HP')
        theta_oc3h = self.ff_reader.get_angle_t0('OP-C3-HP')
        theta_hc3h = self.ff_reader.get_angle_t0('HP-C3-HP')
        plane_c3h = r_c3h*np.sin(0.5*theta_oc3h)
        perp_c3h = r_c3h*np.cos(0.5*theta_oc3h)

        r_c3o = self.ff_reader.get_bond_r0('C3-OP')
        plane_oc3 = r_c3o*np.sin(0.5*theta_oc3h)
        perp_oc3 = r_c3o*np.cos(0.5*theta_oc3h)

        # CH3O
        self.position[0] = np.array([0.0, 0.0, 0.0])
        self.position[3] = np.array([perp_c3h, 0, plane_c3h])
        self.position[4] = self.position[3] + np.array([-perp_oc3, 0, plane_oc3])
        b_val = r_c3h*np.sin(0.5*theta_hc3h)
        a_val = np.sqrt(r_c3h**2-b_val**2)+perp_c3h
        self.position[1] = np.array([a_val, b_val, plane_c3h])
        self.position[2] = np.array([a_val, -b_val, plane_c3h])

        # case repeats = 1
        if self.repeats == 1:
            self.position[-1:-5:-1, 0:2] = self.position[0:4, 0:2]
            self.position[-1:-5:-1, 2] = -self.position[0:4, 2] + 2*self.position[4, 2]
            return

        r_c2c2 = self.ff_reader.get_bond_r0('C2-C2')
        r_c2h = self.ff_reader.get_bond_r0('C2-HP')
        # theta_oc2c2 = self.ff_reader.get_angle_t0('C2-HP')
        theta_hc2h = self.ff_reader.get_angle_t0('HP-C2-HP')
        theta_oc2h = self.ff_reader.get_angle_t0('OP-C2-HP')
        plane_c2c2 = r_c2c2 * np.sin(0.5 * theta_oc2h)
        perp_c2c2 = r_c2c2 * np.cos(0.5 * theta_oc2h)
        plane_c2h = r_c2h * np.sin(0.5 * theta_oc2h)
        perp_c2h = r_c2h * np.cos(0.5 * theta_oc2h)

        r_c2o = self.ff_reader.get_bond_r0('C2-OP')
        plane_oc2 = r_c2o * np.sin(0.5 * theta_oc2h)
        perp_oc2 = r_c2o * np.cos(0.5 * theta_oc2h)

        pos_h = np.zeros([2, 3])
        pos_h_inv = np.zeros([2, 3])
        b_val = r_c2h * np.sin(0.5 * theta_hc2h)
        a_val = np.sqrt(r_c2h ** 2 - b_val ** 2)
        pos_h[0] = np.array([a_val, b_val, 0.0])
        pos_h[1] = np.array([a_val, -b_val, 0.0])
        pos_h_inv[1] = np.array([-a_val, -b_val, 0.0])
        pos_h_inv[0] = np.array([-a_val, b_val, 0.0])
        pos_c2c2 = np.array([perp_c2c2, 0, plane_c2c2])
        pos_c2c2_inv = np.array([-perp_c2c2, 0, plane_c2c2])
        pos_oc2 = np.array([perp_oc2, 0, plane_oc2])
        pos_oc2_inv = np.array([-perp_oc2, 0, plane_oc2])

        for ind_p in range(0, self.repeats-1):
            ind_num = 5+7*ind_p
            if ind_p % 2 == 0:
                self.position[ind_num] = self.position[ind_num-1] + pos_oc2
                self.position[ind_num+1:ind_num+3] = self.position[ind_num] + pos_h[:]
                self.position[ind_num+5] = self.position[ind_num] + pos_c2c2_inv
                self.position[ind_num+3:ind_num+5] = self.position[ind_num+5] + pos_h_inv[:]
                self.position[ind_num+6] = self.position[ind_num+5] + pos_oc2
            else:
                self.position[ind_num] = self.position[ind_num-1] + pos_oc2_inv
                self.position[ind_num+1:ind_num+3] = self.position[ind_num] + pos_h_inv[:]
                self.position[ind_num+5] = self.position[ind_num] + pos_c2c2
                self.position[ind_num+3:ind_num+5] = self.position[ind_num+5] + pos_h[:]
                self.position[ind_num+6] = self.position[ind_num+5] + pos_oc2_inv

        num = 4 + 7 * (self.repeats-1)
        self.position[-1:-5:-1, 2] = -self.position[0:4, 2] + self.position[4, 2] + self.position[num, 2]
        if self.repeats % 2 == 1:
            self.position[-1:-5:-1, 0:2] = self.position[0:4, 0:2]
        else:
            delt_vec = np.array([-0.0184695, 0])
            self.position[-1:-5:-1, 0:2] = -self.position[0:4, 0:2] + self.position[num, 0:2] + delt_vec[:]


    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[-1] - self.position[0]

    def get_bonds(self):
        """
        override of method in SolventAbs
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system
        
        :return: a (self.repeats - 1)x2 numpy array
        """

        b_list = [[3, 4], [5, 4]]
        if self.repeats > 1:
            for ind in range(self.repeats-1):
                b_list += [[5+7*ind, 10+7*ind], [10+7*ind, 11+7*ind], [12+7*ind, 11+7*ind]]

        b_list = [[n+3 for n in sub] for sub in b_list] # changed
        b_list += [[6, 3]] # changed
        b_list += [[3, 0]] #, [3, 1], [3, 2]] # changed

        return np.array(b_list, dtype=int)

    def get_bond_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        
        :return: list of all the bond types in the system
        """
        b_types = ['C2-OP']
        if self.repeats == 1:
            b_types += ['C3-OP']
        else:
            b_types += ['C2-OP']+(self.repeats-2)*['C2-C2', 'C2-OP', 'C2-OP']
            b_types += ['C2-C2'] + ['C2-OP'] + ['C3-OP']

        b_types += ['C2-CS2'] # changed
        b_types += ['CS2-SD'] #+ 2*['CS2-HA'] # changed

        return b_types

    def get_angles(self):
        """
        override of method in SolventAbs
        used to dump_gsd
        
        :return: a list of particle indices (3 particles per angle) involved in the angles in the system
        """

        a_list = [[1, 3, 2], [4, 3, 1], [4, 3, 2], [5, 4, 3]]
        if self.repeats == 1:
            a_list += [[4, 5, 6], [4, 5, 7], [4, 5, 8], [6, 5, 7], [6, 5, 8], [7, 5, 8]]
        else:
            for ind in range(self.repeats-1):
                a_num = 5 + 7*ind
                a_list += [[a_num-1, a_num, a_num+1], [a_num-1, a_num, a_num+2], [a_num+1, a_num, a_num+2]] # [o-c2-h, o-c2-h, h-c2-h]
                a_list += [[a_num-1, a_num, a_num+5], [a_num+3, a_num+5, a_num+4], [a_num+6, a_num+5, a_num]] # [o-c-c, h-c2-h, o-c-c]
                a_list += [[a_num+6, a_num+5, a_num+3], [a_num+6, a_num+5, a_num+4]] # [o-c2-h, o-c2-h]
                a_list += [[a_num, a_num+5, a_num+3], [a_num, a_num+5, a_num+4], [a_num+5, a_num, a_num+1], [a_num+5, a_num, a_num+2]] # [c2-c2-h]*4
                a_list += [[a_num+5, a_num+6, a_num+7]] # [c2-o-c2]
            a_list = a_list[:-1]
            a_list += [[a_num+7, a_num+6, a_num+5]] # [c3-o-c2]
            a_num += 7
            a_list += [[a_num-1, a_num, a_num+1], [a_num-1, a_num, a_num+2], [a_num-1, a_num, a_num+3]]
            a_list += [[a_num+1, a_num, a_num+2], [a_num+2, a_num, a_num+3], [a_num+1, a_num, a_num+3]]

        a_list = [[n+3 for n in sub] for sub in a_list]
        a_list += [[3, 6, 4], [3, 6, 5]] + [[3, 6, 7]]
        a_list += [[0, 3, 1], [0, 3, 2]] + [[1, 3, 2]]
        a_list += [[6, 3, 1], [6, 3, 2]] + [[0, 3, 6]]

        return np.array(a_list, dtype=int)

    def get_angle_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        
        :return: list of all the angle types in the system
        """

        a_types = 1*['HP-C2-HP'] + 2*['OP-C2-HP'] + ['C3-OP-C2'] # changed
        if self.repeats == 1:
            a_types += 3*['OP-C3-HP'] + 3*['HP-C3-HP']
        else:
            l_types = 2*['OP-C2-HP'] + ['HP-C2-HP'] + ['OP-C2-C2'] + ['HP-C2-HP'] + ['OP-C2-C2']
            l_types += 2*['OP-C2-HP']
            l_types += 4*['C2-C2-HP']
            l_types += ['C2-OP-C2']
            a_types += (self.repeats-1)*l_types
            a_types = a_types[:-1]
            a_types += ['C3-OP-C2'] + 3*['OP-C3-HP'] + 3*['HP-C3-HP']

        a_types += 2*['CS2-C2-HP'] + ['CS2-C2-OP']
        a_types += 2*['SD-CS2-HA'] + ['HA-CS2-HA']
        a_types += 2*['C2-CS2-HA'] + ['SD-CS2-C2']

        return a_types

    def get_dihedrals(self):
        """
        override of method in SolventAbs
        used to dump_gsd
        
        :return: a list of particle indices (4 particles per dihedral) involved in the dihedrals in the system
        """

        if self.repeats == 1:
            d_list = [[5, 4, 3, 2], [5, 4, 3, 1], [0, 3, 4, 5]] # changed
            d_list += [[3, 4, 5, 6], [3, 4, 5, 7], [3, 4, 5, 8]]
        else:
            d_list = [[5, 4, 3, 2], [5, 4, 3, 1], [0, 3, 4, 5]] # changed
            d_list += [[3, 4, 5, 6], [3, 4, 5, 7], [3, 4, 5, 10]] # [c3-o-c2-c2]
            d_num = 5
            for ind in range(self.repeats-2):
                d_num = 5 + 7*(ind+1)
                d_list += [[d_num-2, d_num-1, d_num, d_num+1], [d_num-2, d_num-1, d_num, d_num+2]] # [c-o-c-h]
                d_list += [[d_num, d_num-1, d_num-2, d_num-4], [d_num, d_num-1, d_num-2, d_num-3]] # [c-o-c-h]
                d_list += [[d_num-2, d_num-1, d_num, d_num+5], [d_num, d_num-1, d_num-2, d_num-7]] # [c-o-c-c]
            d_num += 7
            d_list += [[d_num, d_num-1, d_num-2, d_num-7]] # [c3-o-c2-c2]
            d_list += [[d_num, d_num-1, d_num-2, d_num-3], [d_num, d_num-1, d_num-2, d_num-4]] # [c3-o-c2-h]
            d_list += [[d_num-2, d_num-1, d_num, d_num+1], [d_num-2, d_num-1, d_num, d_num+2], [d_num-2, d_num-1, d_num, d_num+3]] # [c2-o-c3-h]

            for ind in range(self.repeats-1):
                d_num = 5 + 7*ind
                d_list += [[d_num+6, d_num+5, d_num, d_num+1], [d_num+6, d_num+5, d_num, d_num+2]] # [o-c-c-h]
                d_list += [[d_num-1, d_num, d_num+5, d_num+3], [d_num-1, d_num, d_num+5, d_num+4]] # [o-c-c-h]
                d_list += [[d_num+1, d_num, d_num+5, d_num+3], [d_num+2, d_num, d_num+5, d_num+3]] # [h-c-c-h]
                d_list += [[d_num+1, d_num, d_num+5, d_num+4], [d_num+2, d_num, d_num+5, d_num+4]] # [h-c-c-h]
                d_list += [[d_num-1, d_num, d_num+5, d_num+6]] # [o-c-c-o]

        d_list = [[n+3 for n in sub] for sub in d_list]
        d_list += [[0, 3, 6, 7]] + [[7, 6, 3, 1], [7, 6, 3, 2]]
        d_list += [[4, 6, 3, 0], [5, 6, 3, 0]] + [[4, 6, 3, 1], [4, 6, 3, 2], [5, 6, 3, 1], [5, 6, 3, 2]]


        return np.array(d_list, dtype=int)

    def get_dihedral_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd

        :return: list of all the dihedral types in the system
        """

        d_types = []
        if self.repeats == 1:
            d_types += 2*['C3-OP-C2-HP'] + ['CS2-C2-OP-C3'] + 3*['C2-OP-C3-HP'] # changed
        else:
            d_types += 2*['C2-OP-C2-HP'] + ['CS2-C2-OP-C2'] + 2*['C2-OP-C2-HP'] + ['C2-OP-C2-C2']  # changed
            l_types = 4*['C2-OP-C2-HP']
            l_types += 2*['C2-OP-C2-C2']
            d_types += (self.repeats-2)*l_types
            d_types += ['C3-OP-C2-C2'] + 2*['C3-OP-C2-HP'] + 3*['C2-OP-C3-HP']
            l_types = 4*['OP-C2-C2-HP']
            l_types += 4*['HP-C2-C2-HP']
            l_types += ['OP-C2-C2-OP']
            d_types += (self.repeats-1)*l_types


        d_types += ['SD-CS2-C2-OP'] + 2*['OP-C2-CS2-HA']
        d_types += 2*['HP-C2-CS2-SD'] + 4*['HP-C2-CS2-HA']

        return d_types

    def get_pairs(self):
        """
        override of method in SolventAbs
        used to dump_gsd
        
        :return: a array of 1-4 pairs (particle indices of first and fourth particle for all dihedrals)
        """

        s_list = []
        for d in self.get_dihedrals():
            s_list.append([d[0], d[3]])

        return np.array(s_list, dtype=int)

    def get_pair_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd

        :return: list of all the dihedral types in the system
        """

        return self.get_dihedral_types()

    def get_constraints(self):
        """
        override of method in SolventAbs
        returns a list of the particle indices (2 particles per constraint) involved in the constraints in the system
        
        :return: a (self.repeats - 1)x2 numpy array
        """

        c_list = [[3, 1], [3, 2]]
        if self.repeats == 1:
            c_list += [[5, 6], [5, 7], [5, 8]]
        else:
            for ind in range(self.repeats-1):
                num = 5 + 7*ind
                c_list += [[num, num+1], [num, num+2], [num+5, num+3], [num+5, num+4]]
            num = 5+7*(self.repeats-1)
            c_list += [[num, num+1], [num, num+2], [num, num+3]]

        c_list = [[n+3 for n in sub] for sub in c_list] # changed, added
        c_list += [[3, 1], [3, 2]] # changed, added

        return np.array(c_list, dtype=int)

    def get_constraint_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        
        :return: list of all the constraint types in the system
        """

        c_types = 2*['C2-HP'] # changed
        if self.repeats > 1:
            c_types += 4*(self.repeats-1)*['C2-HP']
        c_types += 3*['C3-HP']
        c_types += 2*['CS2-HA'] # changed, added

        dist_c = []
        for ctyp in c_types:
            dist_c.append(self.ff_reader.get_bond_r0(ctyp))

        return dist_c
