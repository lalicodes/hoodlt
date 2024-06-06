"""
:module: OleateLigand
:platform: Unix, Windows
:synopsis: Implements the oleate (a salt or ester of oleic acid) version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> March 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
..                Xun Zha <xzha@iastate.edu> Sept. 2020
..                  - minor corrections
..                Alex Travesset <trvsst@ameslab.gov> June 2021
                    - Fixed documentation issues
"""

import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class OleateLigand(LigandAbs):
    """
    Defines the oleate (a salt or ester of oleic acid) ligand.
    """

    def __init__(self, forcefield, straight=False):
        """

        :param forcefield: name of the forcefield
        :param straight: if true, initial positions of the ligand form a straight line
        """

        super(OleateLigand, self).__init__(16, forcefield, 20, 'Oleate')

        # particle types
        self.types = ['Coo', 'O', 'CH2', 'CH', 'CH3']
        self.typeid = ['Coo'] + ['O'] * 2 + ['CH2'] * 7 + ['CH'] * 2 + ['CH2'] * 7 + ['CH3']

        self.mass = [self.ff_reader.get_molecular_weight('Coo')] + \
                    [self.ff_reader.get_molecular_weight('O')]*2 + \
                    [self.ff_reader.get_molecular_weight('CH2')]*7 + \
                    [self.ff_reader.get_molecular_weight('CH')]*2 + \
                    [self.ff_reader.get_molecular_weight('CH2')]*7 + \
                    [self.ff_reader.get_molecular_weight('CH3')]
        self.mass = np.array(self.mass)

        # particle positions (Angstroms)
        self.position = self.build_chain(straight)

    def build_chain(self, straight):
        """
        called by initializer to get the positions

        :param straight: if true, initial positions of the ligand form a straight line
        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        o_bond_length = self.ff_reader.get_bond_r0('Coo-O')
        g_bond_length = self.ff_reader.get_bond_r0('Coo-CH2')
        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        m_bond_length = self.ff_reader.get_bond_r0('CH2-CH')
        d_bond_length = self.ff_reader.get_bond_r0('CH-CH')
        o_angle = self.ff_reader.get_angle_t0('O-Coo-O')
        g_angle = self.ff_reader.get_angle_t0('Coo-CH2-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')
        m_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH')
        d_angle = self.ff_reader.get_angle_t0('CH2-CH-CH')

        chain[1] = np.array([-np.cos(o_angle)/2., 0.0, np.sin(o_angle)/2.])*o_bond_length
        chain[2] = np.array([-np.cos(o_angle)/2., 0.0, -np.sin(o_angle)/2.])*o_bond_length
        chain[3:] += np.array([np.sin(g_angle-c_angle/2.), np.cos(g_angle-c_angle/2.), 0.0])*g_bond_length

        ch2vec1 = np.array([np.sin(c_angle/2.0), np.cos(c_angle/2.0), 0.0])*c_bond_length
        ch2vec2 = np.array([np.sin(c_angle/2.0), -np.cos(c_angle/2.0), 0.0])*c_bond_length

        if straight:
            for i in range(4, len(chain)):
                if i % 2 == 0:
                    chain[i:] += ch2vec2
                else:
                    chain[i:] += ch2vec1
        else:
            for i in range(4, 10):
                if i % 2 == 0:
                    chain[i:] += ch2vec2
                else:
                    chain[i:] += ch2vec1

            chain[10:] += np.array([np.sin(m_angle-c_angle/2.), -np.cos(m_angle-c_angle/2.), 0.0])*m_bond_length
            t_ang = d_angle-m_angle+c_angle/2.
            chain[11:] += np.array([np.sin(t_ang), np.cos(t_ang), 0.0])*d_bond_length
            t_ang = d_angle-np.pi+t_ang
            chain[12:] += np.array([np.sin(t_ang), np.cos(t_ang), 0.0])*m_bond_length

            t_ang = m_angle-t_ang
            ch2vec1p = np.array([np.sin(t_ang), -np.cos(t_ang), 0.0])*c_bond_length
            ch2vec2p = np.array([np.sin(c_angle-t_ang), np.cos(c_angle-t_ang), 0.0])*c_bond_length
            for i in range(13, len(chain)):
                if i % 2 == 0:
                    chain[i:] += ch2vec2p
                else:
                    chain[i:] += ch2vec1p

        return chain

    @staticmethod
    def get_bond_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['Coo-O']*2 + ['Coo-CH2'] + ['CH2-CH2']*6 + ['CH2-CH'] + ['CH-CH'] + ['CH2-CH'] + ['CH2-CH2']*7

    @staticmethod
    def get_angle_types():
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        angle = ['O-Coo-O'] + ['O-Coo-CH2']*2 + ['Coo-CH2-CH2'] + ['CH2-CH2-CH2']*5 + ['CH2-CH2-CH'] + \
                ['CH2-CH-CH']*2 + ['CH2-CH2-CH'] + ['CH2-CH2-CH2']*6
        return angle

    @staticmethod
    def get_dihedral_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        dihedral = ['O-Coo-CH2-CH2']*2 + ['Coo-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2']*4 + ['CH2-CH2-CH2-CH'] + \
                   ['CH2-CH2-CH-CH'] + ['CH2-CH-CH-CH2'] + ['CH2-CH2-CH-CH'] + ['CH2-CH2-CH2-CH'] + \
                   ['CH2-CH2-CH2-CH2']*5
        return dihedral

    def get_bonds(self):
        """
        override of the method in LigandAbs
        used by other classes to add bonds when the gsd is dumped
        :return: an n by 2 list of index pairs that are bonded together
        """

        bonds = np.zeros((self.num_particles-1, 2))
        bonds[:3] = np.array([[0, 1], [0, 2], [0, 3]])
        for i in range(3, self.num_particles-1):
            bonds[i] = [i, i+1]
        return bonds

    def get_angles(self):
        """
        override of the method in LigandAbs
        used by other classes to add angles when the gsd is dumped
        :return: an n by 3 list of indexes that have an angle together
        """
        angles = np.zeros((self.num_particles-1, 3))
        angles[:4] = np.array([[1, 0, 2], [1, 0, 3], [2, 0, 3], [0, 3, 4]])
        for i in range(3, self.num_particles-2):
            angles[i+1] = [i, i+1, i+2]
        return angles

    def get_dihedrals(self):
        """
        override of the method in LigandAbs
        used by other classes to add dihedrals when the gsd is dumped
        :return: an n by 4 list of indexes that have a dihedral together
        """
        dis = np.zeros((self.num_particles-3, 4))
        dis[:3] = np.array([[1, 0, 3, 4], [2, 0, 3, 4], [0, 3, 4, 5]])
        for i in range(3, self.num_particles-3):
            dis[i] = [i, i+1, i+2, i+3]
        return dis
