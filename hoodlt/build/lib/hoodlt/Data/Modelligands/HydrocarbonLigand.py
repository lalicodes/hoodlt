"""
:module: HydorcarbonLigand
:platform: Unix, Windows
:synopsis: Implements the hydrocarbon version of the abstract ligand class 

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> April 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class HydrocarbonLigand(LigandAbs):
    """
    Defines the Hydrocarbon ligand
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: number of CH2 in the chain. Total length will have two more, S and CH3
        :param forcefield: string which tells which forcefield to build the object from
        """

        super(HydrocarbonLigand, self).__init__(repeats, ff, repeats+2, 'Hydrocarbon')

        self.mass[0] = self.ff_reader.get_molecular_weight('S')
        self.mass[1:repeats + 1] = self.ff_reader.get_molecular_weight('CH2')
        self.mass[repeats + 1] = self.ff_reader.get_molecular_weight('CH3')

        # types (CH2, CH3)
        self.types = ['S', 'CH2', 'CH3']
        self.typeid = ['S']*1 + ['CH2']*repeats + ['CH3']*1

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions
        
        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        c_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        s_bond_length = self.ff_reader.get_bond_r0('S-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')
        s_angle = self.ff_reader.get_angle_t0('S-CH2-CH2')

        ch2vec1 = np.array([c_bond_length*np.sin(c_angle/2.0), c_bond_length*np.cos(c_angle/2.0), 0.0])
        ch2vec2 = np.array([c_bond_length*np.sin(c_angle/2.0), -c_bond_length*np.cos(c_angle/2.0), 0.0])

        s_pos = np.array([-s_bond_length*np.sin(s_angle/2.0), s_bond_length*np.cos(s_angle/2.0), 0.0])

        current = np.array([0.0, 0.0, 0.0])
        chain[0] = s_pos
        chain[1] = current

        for i in range(2, self.repeats + 2):
            if i % 2 == 0:
                current = np.add(current, ch2vec1)
            else:
                current = np.add(current, ch2vec2)

            chain[i] = current

        for i in range(self.repeats + 2):
            chain[i] = np.subtract(chain[i], s_pos)

        return chain

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        
        :return: list of all the angle types in the chain 
        """
        return ['S-CH2-CH2'] + ['CH2-CH2-CH2'] * (self.repeats - 1)

    def get_bond_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system 
        """
        return ['S-CH2'] + ['CH2-CH2'] * self.repeats

    def get_dihedral_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """
    
        # there are two returns here because the spec one is hack
        # return ['CH2-CH2-CH2-CH2-spec'] * (self.repeats - 1)
        return ['CH2-CH2-CH2-CH2'] * (self.repeats - 1)
