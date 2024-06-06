"""
:module: Alkane
:platform: Unix, Windows
:synopsis: Implements a class defining Alkane solvents

.. moduleauthor:: John Mobley IV <jmobley@iastate.edu> June 2018
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


class Alkane(SolventAbs):
    """
    Defines all Alkane solvents
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        super(Alkane, self).__init__(repeats, ff, repeats, 'Alkane')

        # in the future, we would like to be able to build methane
        if repeats < 2:
            raise ValueError("Number of repeats must be greater than or equal to 2")

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('CH3')
        self.mass[1:repeats - 1] = self.ff_reader.get_molecular_weight('CH2')
        self.mass[repeats - 1] = self.ff_reader.get_molecular_weight('CH3')

        # types (CH2, CH3)
        self.types = ['CH2', 'CH3']
        self.typeid = ['CH3']*1 + ['CH2']*(repeats - 2) + ['CH3']*1

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions
        
        :return: the position array for the Alkane chain
        """

        chain = np.zeros((self.num_particles, 3))

        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')

        ch2vec1 = np.array([c_bond_length*np.sin(c_angle/2.0), c_bond_length*np.cos(c_angle/2.0), 0.0])
        ch2vec2 = np.array([c_bond_length*np.sin(c_angle/2.0), -c_bond_length*np.cos(c_angle/2.0), 0.0])

        ch3_pos = np.array([-c_bond_length*np.sin(c_angle/2.0), c_bond_length*np.cos(c_angle/2.0), 0.0])

        current = np.array([0.0, 0.0, 0.0])
        chain[0] = ch3_pos
        chain[1] = current

        for i in range(2, self.repeats):
            if i % 2 == 0:
                current = np.add(current, ch2vec1)
            else:
                current = np.add(current, ch2vec2)

            chain[i] = current

        for i in range(self.repeats):
            chain[i] = np.subtract(chain[i], ch3_pos)

        return chain

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

        b_list = []
        i = 0
        while i < self.repeats:
            if i:
                b_list.append([i - 1, i])
            i = i + 1

        return np.array(b_list)

    def get_angles(self):
        """
        override of method in SolventAbs
        returns a list of particle indices (3 particles per angle) involved in the angles in the system
        
        :return: a (self.repeats - 2)x3 numpy array
        """

        a_list = []
        i = 0
        while i < self.repeats - 1:
            if i:
                a_list.append([i - 1, i, i + 1])
            i = i + 1
        return np.array(a_list)

    def get_dihedrals(self):
        """
        override of method in SolventAbs
        returns a list of particle indices (4 particles per dihedral) involved in the dihedrals in the system
        
        :return: a (self.repeats - 3)x4 numpy array
        """

        d_list = []
        i = 0
        while i < self.repeats - 2:
            if i:
                d_list.append([i - 1, i, i + 1, i + 2])
            i = i + 1
        return np.array(d_list)

    def get_bond_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        
        :return: list of all the bond types in the system
        """

        return ['CH2-CH2'] * (self.repeats - 1)

    def get_angle_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd
        
        :return: list of all the angle types in the system
        """

        return ['CH2-CH2-CH2'] * (self.repeats - 2)

    def get_dihedral_types(self):
        """
        override of the method in SolventAbs
        used to dump_gsd

        :return: list of all the dihedral types in the system
        """

        return ['CH2-CH2-CH2-CH2'] * (self.repeats - 3)
