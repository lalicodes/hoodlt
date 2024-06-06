"""
:module: OleylamineLigand
:platform: Unix, Windows
:synopsis: Implements the oleyl ammonium version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> March 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
..                Xun Zha <xzha@iastate.edu> Sept. 2020
..                  - minor corrections
..                Xun Zha <xzha@iastate.edu> Dec. 2020
..                  - corrections
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class OleylammoniumLigand(LigandAbs):
    """
    Defines the oleyl ammonium ligand
    """

    def __init__(self, forcefield, straight=False):
        """

        :param forcefield: name of the forcefield
        :param straight: if true, initial positions of the ligand form a straight line
        """

        super(OleylammoniumLigand, self).__init__(17, forcefield, 19, 'Oleylammonium')

        # particle types
        self.types = ['NH3', 'CH2', 'CH', 'CH3']
        self.typeid = ['NH3'] + ['CH2'] * 8 + ['CH'] * 2 + ['CH2'] * 7 + ['CH3']

        # get units and FF parameters

        # molecular weight and charge
        self.mass = [self.ff_reader.get_molecular_weight('NH3')] + \
                    [self.ff_reader.get_molecular_weight('CH2')] * 8 + \
                    [self.ff_reader.get_molecular_weight('CH')] * 2 + \
                    [self.ff_reader.get_molecular_weight('CH2')] * 7 + \
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

        n_bond_length = self.ff_reader.get_bond_r0('NH3-CH2')
        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        m_bond_length = self.ff_reader.get_bond_r0('CH2-CH')
        d_bond_length = self.ff_reader.get_bond_r0('CH-CH')
        n_angle = self.ff_reader.get_angle_t0('NH3-CH2-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')
        m_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH')
        d_angle = self.ff_reader.get_angle_t0('CH2-CH-CH')

        chain[1:] += np.array([np.sin(n_angle/2.), np.cos(n_angle/2.), 0.0])*n_bond_length

        ch2vec1 = np.array([np.sin(c_angle/2.0), np.cos(c_angle/2.0), 0.0])*c_bond_length
        ch2vec2 = np.array([np.sin(c_angle/2.0), -np.cos(c_angle/2.0), 0.0])*c_bond_length

        if straight:
            for i in range(2, len(chain)):
                if i % 2 == 0:
                    chain[i:] += ch2vec2
                else:
                    chain[i:] += ch2vec1
        else:
            for i in range(2, 9):
                if i % 2 == 0:
                    chain[i:] += ch2vec2
                else:
                    chain[i:] += ch2vec1

            chain[9:] += np.array([np.sin(m_angle-c_angle/2.), np.cos(m_angle-c_angle/2.), 0.0])*m_bond_length
            t_ang = d_angle-m_angle+c_angle/2.
            chain[10:] += np.array([np.sin(t_ang), -np.cos(t_ang), 0.0])*d_bond_length
            t_ang = d_angle-np.pi+t_ang
            chain[11:] += np.array([np.sin(t_ang), -np.cos(t_ang), 0.0])*m_bond_length

            t_ang = m_angle-t_ang
            ch2vec1p = np.array([np.sin(t_ang), np.cos(t_ang), 0.0])*c_bond_length
            ch2vec2p = np.array([np.sin(c_angle-t_ang), -np.cos(c_angle-t_ang), 0.0])*c_bond_length
            for i in range(12, len(chain)):
                if i % 2 == 0:
                    chain[i:] += ch2vec1p
                else:
                    chain[i:] += ch2vec2p

        return chain

    @staticmethod
    def get_bond_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['NH3-CH2'] + ['CH2-CH2']*7 + ['CH2-CH'] + ['CH-CH'] + ['CH2-CH'] + ['CH2-CH2']*7

    @staticmethod
    def get_angle_types():
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        angle = ['NH3-CH2-CH2'] + ['CH2-CH2-CH2']*6 + ['CH2-CH2-CH'] + ['CH2-CH-CH']*2 + ['CH2-CH2-CH'] + \
                ['CH2-CH2-CH2']*6
        return angle

    @staticmethod
    def get_dihedral_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        dihedral = ['NH3-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2']*5 + ['CH2-CH2-CH2-CH'] + ['CH2-CH2-CH-CH'] + \
                   ['CH2-CH-CH-CH2'] + ['CH2-CH2-CH-CH'] + ['CH2-CH2-CH2-CH'] + ['CH2-CH2-CH2-CH2']*5
        return dihedral
