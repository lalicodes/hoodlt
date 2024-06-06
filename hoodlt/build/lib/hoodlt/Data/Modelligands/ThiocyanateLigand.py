"""
:module: ThiocyanateLigand
:platform: Unix, Windows
:synopsis: Implements the thiocyanate (rhodanide) version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> December 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
..                Xun Zha <xzha@iastate.edu> Sept. 2020
..                  - minor corrections
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs

class ThiocyanateLigand(LigandAbs):
    """
    Defines the thiocyanate ligand
    """

    def __init__(self, forcefield):
        """

        :param forcefield: name of the forcefield
        """

        super(ThiocyanateLigand, self).__init__(1, forcefield, 3, 'Thiocyanate')  # always 1 repeat

        # types
        self.types = ['S', 'C', 'N']
        self.typeid = ['S', 'C', 'N']

        self.mass = np.array([self.ff_reader.get_molecular_weight('S'),
                              self.ff_reader.get_molecular_weight('C'),
                              self.ff_reader.get_molecular_weight('N')])

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        s_bond_length = self.ff_reader.get_bond_r0('S-C')
        c_bond_length = self.ff_reader.get_bond_r0('C-N')

        chain[1:, 0] += s_bond_length
        chain[2, 0] += c_bond_length

        return chain

    @staticmethod
    def get_bond_types():
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['S-C'] + ['C-N']

    @staticmethod
    def get_angle_types():
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        return ['S-C-N']

    @staticmethod
    def get_dihedrals():
        return []
