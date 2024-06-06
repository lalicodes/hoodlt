"""
:module: AmmoniumLigand
:platform: Unix, Windows
:synopsis: Implements the ammonium version of the abstract ligand class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> January 2018
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


class AmmoniumLigand(LigandAbs):
    """
    Defines the empty ligand
    """

    def __init__(self, forcefield):
        """Defines the ammonium ligand

        :param forcefield: name of the forcefield
        """

        graft_type = 'NH4'
        super(AmmoniumLigand, self).__init__(1, forcefield, 1, graft_type)

        # types
        self.types = [graft_type]
        self.typeid = [graft_type]

        # get values from forcefield
        self.mass = np.array([self.ff_reader.get_molecular_weight(graft_type)])

        # particle positions (Angstroms)
        self.position = np.zeros([1, 3])

    @staticmethod
    def chain_vector():
        """a unit vector
        """
        return np.array([1.0, 0.0, 0.0])

    @staticmethod
    def get_angles():
        return []

    @staticmethod
    def get_dihedrals():
        return []
