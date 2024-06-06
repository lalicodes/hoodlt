"""
:module: CGAtacticPSLigand
:platform: Unix, Windows
:synopsis: Implements the coarse-grained atactic polystyrene version of the abstract ligand class 

.. moduleauthor:: Jianshe Xia <xiajs6075@iccas.ac.cn> December 2018
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class CGAtacticPSLigand(LigandAbs):
    """
    Defines the CGAtacticPS ligand
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: number of aPS monomer in the chain. Total length will have two more, Gr, STH
        :param forcefield: string which tells which forcefield to build the object from
        """

        super(CGAtacticPSLigand, self).__init__(repeats, ff, repeats + 2, 'CGAtacticPS')

        self.mass[0] = self.ff_reader.get_molecular_weight('Gr')
        self.mass[1:repeats + 1] = self.ff_reader.get_molecular_weight('ST')
        self.mass[repeats + 1] = self.ff_reader.get_molecular_weight('STH')

        # types (ST, STH)
        self.types = ['Gr', 'ST', 'STH']
        self.typeid = ['Gr'] * 1 + ['ST'] * repeats + ['STH'] * 1

        # particle positions (Angstroms)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        st_bond_length = self.ff_reader.get_bond_r0('ST-ST')
        gr_bond_length = self.ff_reader.get_bond_r0('Gr-ST')
        st_angle = self.ff_reader.get_angle_t0('ST-ST-ST')
        gr_angle = self.ff_reader.get_angle_t0('Gr-ST-ST')

        ch2vec1 = np.array([st_bond_length * np.sin(st_angle / 2.0), st_bond_length * np.cos(st_angle / 2.0), 0.0])
        ch2vec2 = np.array([st_bond_length * np.sin(st_angle / 2.0), -st_bond_length * np.cos(st_angle / 2.0), 0.0])

        gr_pos = np.array([-gr_bond_length * np.sin(gr_angle / 2.0), gr_bond_length * np.cos(gr_angle / 2.0), 0.0])

        current = np.array([0.0, 0.0, 0.0])
        chain[0] = gr_pos
        chain[1] = current

        for i in range(2, self.repeats + 2):
            if i % 2 == 0:
                current = np.add(current, ch2vec1)
            else:
                current = np.add(current, ch2vec2)

            chain[i] = current

        for i in range(self.repeats + 2):
            chain[i] = np.subtract(chain[i], gr_pos)

        return chain

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain
        """
        return ['Gr-ST-ST'] + ['ST-ST-ST'] * (self.repeats - 1)

    def get_bond_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system
        """
        return ['Gr-ST'] + ['ST-ST'] * self.repeats

    def get_dihedral_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the dihedral types in the system
        """

        # there are two returns here because the spec one is hack
        return ['ST-ST-ST-ST'] * (self.repeats - 1)
