"""
:module: AlkylAmmoniumLigand
:platform: Unix, Windows
:synopsis: Implements the alkyl ammonium version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> December 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class AlkylAmmoniumLigand(LigandAbs):
    """
    Defines the alkyl ligand
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: the number of repeat units in the chain
        """

        super(AlkylAmmoniumLigand, self).__init__(repeats, ff, repeats+2, 'AlkylAmmonium')

        # particle types
        self.types = ['NH3', 'CH2', 'CH3']
        self.typeid = ['NH3'] + ['CH2'] * repeats + ['CH3']

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('NH3')
        self.mass[1:-1] = self.ff_reader.get_molecular_weight('CH2')
        self.mass[-1] = self.ff_reader.get_molecular_weight('CH3')

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        n_bond_length = self.ff_reader.get_bond_r0('NH3-CH2')
        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        n_angle = self.ff_reader.get_angle_t0('NH3-CH2-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')

        chain[1:] += np.array([np.sin(n_angle-c_angle/2.), np.cos(n_angle-c_angle/2.), 0.0])*n_bond_length

        ch2vec1 = np.array([np.sin(c_angle/2.0), np.cos(c_angle/2.0), 0.0])*c_bond_length
        ch2vec2 = np.array([np.sin(c_angle/2.0), -np.cos(c_angle/2.0), 0.0])*c_bond_length
        for i in range(2, len(chain)):
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
        return ['NH3-CH2'] + ['CH2-CH2'] * self.repeats

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        return ['NH3-CH2-CH2'] + ['CH2-CH2-CH2'] * (self.repeats-1)

    def get_dihedral_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        return ['NH3-CH2-CH2-CH2'] + ['CH2-CH2-CH2-CH2'] * (self.repeats-2)
