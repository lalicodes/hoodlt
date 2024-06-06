"""
:module: CysteineCrystal
:platform: Unix, Windows
:synopsis: Implements a class defining Ethane

.. moduleauthor:: Jacob Austin <jaustin2@iastate.edu> November 2019
.. history:
"""


from __future__ import division
import numpy as np
import importlib_resources
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs

file_with_positions = 'OrthorhombicCysteine.txt'
d_name = 'Data/Modelsolvent/' + file_with_positions
ref = importlib_resources.files('hoodlt')/d_name
with importlib_resources.as_file(ref) as path:
    pos_list = (np.genfromtxt(path)).tolist()

class OrthorhombicCysteine(SolventAbs):
    """
    Defines solvent OrthorhombicCysteine in the orthorhombic configuration
    """

    def __init__(self, ff):
        super(OrthorhombicCysteine, self).__init__(1, ff, 56, 'OrthorhombicCysteine')

        self.types = ['N', 'H_1', 'H_2', 'H_3', 'C_4', 'H_5', 'S', 'C_6', 'C_8', 'O_9', 'O_A', 'H_B', 'H_C', 'H_d']
        self.typeid = ['N', 'H_1', 'H_2', 'H_3', 'C_4', 'H_5', 'S', 'C_6', 'C_8', 'O_9', 'O_A', 'H_B', 'H_C', 'H_d']*4

        # masses
        for i in range(3):
            self.mass[0 + i*14] = self.ff_reader.get_molecular_weight('N')
            self.mass[1 + i*14] = self.ff_reader.get_molecular_weight('H_1')
            self.mass[2 + i*14] = self.ff_reader.get_molecular_weight('H_2')
            self.mass[3 + i*14] = self.ff_reader.get_molecular_weight('H_3')
            self.mass[4 + i*14] = self.ff_reader.get_molecular_weight('C_4')
            self.mass[5 + i*14] = self.ff_reader.get_molecular_weight('H_5')
            self.mass[7 + i*14] = self.ff_reader.get_molecular_weight('C_6')
            self.mass[6 + i*14] = self.ff_reader.get_molecular_weight('S')
            self.mass[8 + i*14] = self.ff_reader.get_molecular_weight('C_8')
            self.mass[9 + i*14] = self.ff_reader.get_molecular_weight('O_9')
            self.mass[10 + i*14] = self.ff_reader.get_molecular_weight('O_A')
            self.mass[11 + i*14] = self.ff_reader.get_molecular_weight('H_B')
            self.mass[12 + i*14] = self.ff_reader.get_molecular_weight('H_C')
            self.mass[13 + i*14] = self.ff_reader.get_molecular_weight('H_d')

        # charges
        for i in range(3):
            self.charge[0 + i*14] = self.ff_reader.get_charge('N')
            self.charge[1 + i*14] = self.ff_reader.get_charge('H_1')
            self.charge[2 + i*14] = self.ff_reader.get_charge('H_2')
            self.charge[3 + i*14] = self.ff_reader.get_charge('H_3')
            self.charge[4 + i*14] = self.ff_reader.get_charge('C_4')
            self.charge[5 + i*14] = self.ff_reader.get_charge('H_5')
            self.charge[7 + i*14] = self.ff_reader.get_charge('C_6')
            self.charge[6 + i*14] = self.ff_reader.get_charge('S')
            self.charge[8 + i*14] = self.ff_reader.get_charge('C_8')
            self.charge[9 + i*14] = self.ff_reader.get_charge('O_9')
            self.charge[10 + i*14] = self.ff_reader.get_charge('O_A')
            self.charge[11 + i*14] = self.ff_reader.get_charge('H_B')
            self.charge[12 + i*14] = self.ff_reader.get_charge('H_C')
            self.charge[13 + i*14] = self.ff_reader.get_charge('H_d')


        # atom positions (nm)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position of lcysteine
        """

        chain = np.array(pos_list)

        chain = chain / 10

        return chain

    def get_vector(self):
        """Computes the length of the chain
        Overrides method in LigandAbs
        :return: vector from first to last position
        :rtype: numpy pos array [x,y,z]
        """
        return np.subtract(self.position[4], self.position[0])

    def get_bonds(self):
        """ overrides method in LigandAbs
        :return:
        """
        init = np.array([[1, 0], [2, 0], [3, 0], [4, 0], [5, 4], [7, 4], [6, 7], [8, 4], [9, 8], [10, 8],
                         [11, 7], [12, 7], [13, 6]])
        mol2 = np.array(init)
        mol3 = np.array(init)
        mol4 = np.array(init)
        for i in range(len(init)):
            mol2[i]+=14
            mol3[i]+=28
            mol4[i]+=42
        bond = np.concatenate((init, mol2, mol3, mol4), axis=0)
        return bond

    def get_bond_types(self):
        """
        :return:
        """
        return ['H_1-N', 'H_2-N', 'H_3-N', 'C_4-N', 'H_5-C_4', 'C_6-C_4', 'S-C_6', 'C_8-C_4', 'O_9-C_8', 'O_A-C_8',
                'H_B-C_6', 'H_C-C_6', 'H_d-S']*4

    def get_angles(self):
        """ overrides method in LigandAbs
        :return: angles
        """
        init = np.array([[1, 0, 2], [1, 0, 3], [1, 0, 4], [0, 4, 5], [0, 4, 7], [4, 7, 6], [0, 4, 8], [4, 8, 9],
                           [4, 8, 10], [4, 7, 11], [4, 7, 12], [7, 6, 13], [6, 7, 11], [2, 0, 3], [3, 0, 4], [2, 0, 4],
                           [7, 4, 8], [6, 7, 12], [5, 4, 7], [9, 8, 10], [5, 4, 7], [11, 7, 12]])
        mol2 = np.array(init)
        mol3 = np.array(init)
        mol4 = np.array(init)
        for i in range(len(init)):
            mol2[i]+=14
            mol3[i]+=28
            mol4[i]+=42
        angles = np.concatenate((init, mol2, mol3, mol4), axis = 0)
        return angles

    def get_angle_types(self):
        """
        :return:
        """
        return ['H_1-N-H_3', 'H_1-N-H_4', 'H_1-N-C_4', 'N-C_4-H_5', 'N-C_4-C_6', 'C_4-C_6-S', 'N-C_4-C_8', 'C_4-C_8-O_9',
                'C_4-O_9-O_A', 'C_4-C_6-H_B', 'C_4-C_6-H_C', 'C_4-S-H_d', 'S-C_6-H_B', 'H_2-N-H_3', 'H_3-N-C_4',
                'H_2-N-C_4', 'C_6-C_4-C_8', 'S-C_6-H_C', 'H_5-C_4-C_6', 'O_9-C_8-O_A', 'H_5-C_4-C_8', 'H_B-C_6-H_C']*4

    def get_dihedrals(self):
        """ overrides method in LigandAbs
        :return: dihedrals
        """

        init = np.array([[8, 4, 7, 6], [8, 4, 0, 3], [8, 4, 0, 2], [8, 4, 0, 1], [7, 4, 0, 2], [7, 4, 0, 3],
                              [7, 4, 0, 1], [12, 7, 4, 8], [11, 7, 4, 8], [11, 7, 4, 5], [12, 7, 4, 5], [11, 7, 4, 0],
                              [12, 7, 4, 0], [5, 4, 0, 3], [5, 4, 0, 1], [5, 4, 0, 2], [13, 6, 7, 4], [13, 6, 7, 11],
                              [13, 6, 7, 12], [6, 7, 4, 5], [6, 7, 4, 0]])
        mol2 = np.array(init)
        mol3 = np.array(init)
        mol4 = np.array(init)
        for i in range(len(init)):
            mol2[i]+=14
            mol3[i]+=28
            mol4[i]+=42
        dihedrals = np.concatenate((init, mol2, mol3, mol4), axis = 0)

        return dihedrals

    def get_dihedral_types(self):
        """

        :return:
        """
        return ['C_8-C_4-C_6-S', 'C_8-C_4-N-H_3', 'C_8-C_4-N-H_2', 'C_8-C_4-N-H_1', 'C_6-C_4-N-H_2', 'C_6-C_4-N-H_3',
                'C_6-C_4-N-H_1', 'H_C-C_6-C_4-C_8', 'H_B-C_6-C_4-C_8', 'H_B-C_6-C_4-H_5', 'H_C-C_6-C_4-H_5',
                'H_B-C_6-C_4-N', 'H_C-C_6-C_4-N', 'H_5-C_4-N-H_3', 'H_5-C_4-N-H_1', 'H_5-C_4-N-H_2', 'H_d-S-C_6-H_4',
                'H_d-S-C_6-H_B', 'H_d-S-C_6-H_C', 'S-C_6-C_4-H_5', 'S-C_6-C_4-N']*4
    def get_impropers(self):
        """ overrides method in LigandAbs
        :return:
        """
        impropers = np.array([[10, 8, 4, 9], [24, 22, 18, 23], [38, 36, 32, 37], [52, 40, 46, 53]])
        return impropers

    def get_improper_types(self):

        return ['O_A-C_8-C_4-O_9']*4


