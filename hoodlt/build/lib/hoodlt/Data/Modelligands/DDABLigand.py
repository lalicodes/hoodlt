"""
:module: DDABLigand
:platform: Unix, Windows
:synopsis: Implements the didodecyldimethylammonium bromide (DDAB) version of the abstract ligand class

.. moduleauthor:: Jonas Hallstrom <jlhallst@asu.edu> June 2021
.. history:
..                Jonas Hallstrom <jlhallst@asu.edu> June 2021
..                  - created class based primarily on the OleateLigand.py class by Xun Zha
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class DDABLigand(LigandAbs):
    """
    Defines the DDAB ligand, based on the oleate ligand.
    """

    def __init__(self, forcefield, straight=False):
        """

        :param forcefield: name of the forcefield
        :param straight: if true, initial positions of the ligand form a straight line
        """

        super(DDABLigand, self).__init__(22, forcefield, 27, 'DDAB')

        # particle types
        self.types = ['N', 'CH3', 'CH2']
        self.typeid = ['N'] + ['CH3'] * 2 + ['CH2'] * 11 + ['CH3'] + ['CH2'] * 11 + ['CH3']

        self.mass = [self.ff_reader.get_molecular_weight('N')] + \
                    [self.ff_reader.get_molecular_weight('CH3')]*2 + \
                    [self.ff_reader.get_molecular_weight('CH2')]*11 + \
                    [self.ff_reader.get_molecular_weight('CH3')] + \
                    [self.ff_reader.get_molecular_weight('CH2')]*11 + \
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

        nch3_bond_length = self.ff_reader.get_bond_r0('N-CH3')
        nch2_bond_length = self.ff_reader.get_bond_r0('N-CH2')
        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        ch3_angle = self.ff_reader.get_angle_t0('CH3-N-CH3')
        ch2_angle = self.ff_reader.get_angle_t0('CH2-N-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')

        chain[1] = np.array([-np.cos(ch3_angle/2.),  0.0, np.sin(ch3_angle/2.)])*nch3_bond_length
        chain[2] = np.array([-np.cos(ch3_angle/2.),  0.0, -np.sin(ch3_angle/2.)])*nch3_bond_length
        chain[3:15] += np.array([np.cos(ch2_angle/2.), -np.sin(ch2_angle/2.), 0.0])*nch2_bond_length 
        chain[15:] += np.array([np.cos(ch2_angle/2.), np.sin(ch2_angle/2.), 0.0])*nch2_bond_length

        ch2vec1 = np.array([np.sin(c_angle/2.0), np.cos(c_angle/2.0), 0.0])*c_bond_length
        ch2vec2 = np.array([np.sin(c_angle/2.0), -np.cos(c_angle/2.0), 0.0])*c_bond_length

        if straight:
            for i in range(4, 15):
                if i % 2 == 0:
                    chain[i:15] += ch2vec2
                else:
                    chain[i:15] += ch2vec1
            for i in range(16, self.num_particles):
                if i % 2 == 0:
                    chain[i:] += ch2vec1
                else:
                    chain[i:] += ch2vec2
        else:
            print("There is no non-straight version of the ligand")
                    
        return chain

    @staticmethod
    def get_bond_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['N-CH3']*2 + ['N-CH2'] + ['CH2-CH2']*11 + ['N-CH2'] + ['CH2-CH2']*11

    @staticmethod
    def get_angle_types():
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        angle = ['CH3-N-CH3'] + ['CH2-N-CH2'] + ['CH3-N-CH2']*4 + ['N-CH2-CH2'] + ['CH2-CH2-CH2']*10 + ['N-CH2-CH2'] + ['CH2-CH2-CH2']*10
        return angle

    @staticmethod
    def get_dihedral_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        dihedral = ['CH3-N-CH2-CH2']*4 + ['CH2-N-CH2-CH2']*2 + ['N-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2']*9 + ['N-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2']*9
        return dihedral

    def get_bonds(self):
        """
        override of the method in LigandAbs
        used by other classes to add bonds when the gsd is dumped
        :return: an n by 2 list of index pairs that are bonded together
        """

        bonds = np.zeros((self.num_particles-1, 2))
        bonds[:3] = np.array([[0, 1], [0, 2], [0, 3]])
        for i in range(3, 14):
            bonds[i] = [i, i+1]
        bonds[14]=[0, 15]
        for i in range(15, 26):
            bonds[i] = [i, i+1]        
        return bonds

    def get_angles(self):
        """
        override of the method in LigandAbs
        used by other classes to add angles when the gsd is dumped
        :return: an n by 3 list of indexes that have an angle together
        """
        angles = np.zeros((self.num_particles+1, 3))
        angles[:2] = np.array([[1, 0, 2], [3, 0, 15]]) 
        angles[2:6] = np.array([[1, 0, 3], [1, 0, 15], [2, 0, 3], [2, 0, 15]])
        angles[6] = [0, 3, 4]
        for i in range(3, 13):
            angles[i+4] = [i, i+1, i+2]
        angles[17] = [0, 15, 16]
        for i in range(15, 25):
            angles[i+3] = [i, i+1, i+2]            
        return angles

    def get_dihedrals(self):
        """
        override of the method in LigandAbs
        used by other classes to add dihedrals when the gsd is dumped
        :return: an n by 4 list of indexes that have a dihedral together
        """
        dis = np.zeros((self.num_particles-1, 4))
        dis[:4] = np.array([[1, 0, 3, 4], [1, 0, 15, 16], [2, 0, 3, 4], [2, 0, 15, 16]])
        dis[4:6] = np.array([[15, 0, 3, 4], [3, 0, 15, 16]])
        dis[6] = [0, 3, 4, 5]
        for i in range(3, 12):
            dis[i+4] = [i, i+1, i+2, i+3]
        dis[16] = [0, 15, 16, 17]
        for i in range(15, 24):
            dis[i+2] = [i, i+1, i+2, i+3]            
        return dis
