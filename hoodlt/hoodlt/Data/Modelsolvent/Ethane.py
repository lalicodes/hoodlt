"""
:module: Ethane
:platform: Unix, Windows
:synopsis: Implements a class defining Ethane

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> November 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Removed unecessary functions
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - edited get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs

class Ethane(SolventAbs):
    """
    Defines solvent Ethane  (:math:`\\mbox{CH}_3\\mbox{CH}_3`) 
    """

    def __init__(self, ff):
        """

        :param ff: the name of the forcefield used to build the Ethane
        """

        super(Ethane, self).__init__(2, ff, 2, 'Ethane')  # ethane will always have 2 repeats

        self.mass = np.zeros((2, 1))
        self.mass[0] = self.ff_reader.get_molecular_weight('Ethane')
        self.mass[1] = self.ff_reader.get_molecular_weight('Ethane')

        # types (CH2, CH3)
        self.types = ['Ethane']
        self.typeid = ['Ethane']*2

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position of ethane
        """

        r0 = self.ff_reader.get_bond_r0('CH2-CH2')

        chain = np.array([[0, 0, 0], [r0, 0, 0]])

        return chain

    def get_vector(self):
        """
        See documentation in SolventAbs
        """

        return np.array([1, 0, 0])

    def get_bonds(self):
        """
        override of method in SolventAbs
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system

        :return: a 1x2 numpy array
        """

        return np.array([[0,1]])

    def get_bond_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        :return: list of all the bond types in the system 
        """
        return ['CH2-CH2']
