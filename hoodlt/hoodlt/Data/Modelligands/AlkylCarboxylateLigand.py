"""
:module: AlkylCarboxylateLigand
:platform: Unix, Windows
:synopsis: Implements the alkyl carboxylate version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> December 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class AlkylCarboxylateLigand(LigandAbs):
    """
    Defines the alkyl ligand
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: the number of repeat units in the chain
        """

        super(AlkylCarboxylateLigand, self).__init__(repeats, ff, repeats+4, 'AlkylCarboxylate')

        # particle types
        self.types = ['Coo', 'O', 'CH2', 'CH3']
        self.typeid = ['Coo'] + ['O'] * 2 + ['CH2'] * repeats + ['CH3']

        self.mass[0] = self.ff_reader.get_molecular_weight('Coo')
        self.mass[1:3] = self.ff_reader.get_molecular_weight('O')
        self.mass[3:-1] = self.ff_reader.get_molecular_weight('CH2')
        self.mass[-1] = self.ff_reader.get_molecular_weight('CH3')

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((self.num_particles, 3))

        o_bond_length = self.ff_reader.get_bond_r0('Coo-O')
        g_bond_length = self.ff_reader.get_bond_r0('Coo-CH2')
        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        o_angle = self.ff_reader.get_angle_t0('O-Coo-O')
        g_angle = self.ff_reader.get_angle_t0('Coo-CH2-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')

        chain[1] = np.array([-np.cos(o_angle)/2., 0.0, np.sin(o_angle)/2.])*o_bond_length
        chain[2] = np.array([-np.cos(o_angle)/2., 0.0, -np.sin(o_angle)/2.])*o_bond_length
        chain[3:] += np.array([np.sin(g_angle-c_angle/2.), np.cos(g_angle-c_angle/2.), 0.0])*g_bond_length

        ch2vec1 = np.array([np.sin(c_angle/2.0), np.cos(c_angle/2.0), 0.0])*c_bond_length
        ch2vec2 = np.array([np.sin(c_angle/2.0), -np.cos(c_angle/2.0), 0.0])*c_bond_length
        for i in range(4, len(chain)):
            if i % 2 == 0:
                chain[i:] += ch2vec2
            else:
                chain[i:] += ch2vec1

        return chain

    def get_bond_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['Coo-O'] * 2 + ['Coo-CH2'] + ['CH2-CH2'] * self.repeats

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        return ['O-Coo-O'] + ['O-Coo-CH2'] * 2 + ['Coo-CH2-CH2'] + ['CH2-CH2-CH2'] * (self.repeats-1)

    def get_dihedral_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        return ['O-Coo-CH2-CH2'] * 2 + ['Coo-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2'] * (self.repeats-2)

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
        dis = np.zeros((self.repeats+1, 4))
        dis[:3] = np.array([[1, 0, 3, 4], [2, 0, 3, 4], [0, 3, 4, 5]])
        for i in range(3, 3 + self.repeats-2):
            dis[i] = [i, i+1, i+2, i+3]
        return dis
