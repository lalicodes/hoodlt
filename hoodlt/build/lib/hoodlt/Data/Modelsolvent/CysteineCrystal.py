"""
:module: CysteineCrystal
:platform: Unix, Windows
:synopsis: Implements a class defining Ethane

.. moduleauthor:: Jacob Austin <jaustin2@iastate.edu> November 2019
.. history:
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs

class CysteineCrystal(SolventAbs):
    """
    Defines solvent CysteineCrystal in the monoclinic configuration
    """

    def __init__(self, ff):
        super(CysteineCrystal, self).__init__(1, ff, 56, 'CysteineCrystal')

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

        chain = np.array([[2.20112, -3.95865, -1.33904], [1.53511, -3.61035, -1.88450], [3.06739, -3.87714, -1.71466],
                          [2.16573, -4.94966, -1.30038], [2.06144, -3.26834, -0.03264], [2.75600, -3.59194, 0.54811],
                          [0.83552, -0.88660, -0.94153], [2.23907, -1.75820, -0.18090], [0.71372, -3.62055, 0.60253],
                          [-0.08258, -4.32689, -0.06014], [0.49426, -3.16266,  1.75071], [3.07554, -1.58325, -0.74610],
                          [2.40507, -1.35761,  0.74760], [0.89683, -1.45459, -2.15535], [0.04107, -7.00688, -0.69045],
                          [0.70275, -7.18438, -1.37481], [-0.68538, -7.35225, -0.92747], [-0.15151, -6.05268, -0.62809],
                          [0.46717, -7.53251, 0.63571], [-0.22340, -7.40637, 1.25631], [1.01334, -9.86890, 2.07330],
                          [0.75679, -9.02771, 0.49477], [1.71106, -6.76004, 1.09531], [2.51489, -6.42160, 0.17861],
                          [1.85690, -6.54870, 2.31162], [1.57460, -9.14720, -0.07218], [-0.01002, -9.46222, 0.01746], [-0.32531, -9.94887, 2.26515],
                          [4.94914, -5.88159, 1.22232], [ 5.60031, -6.22232, 0.65451], [4.07308, -5.95815, 0.86902], [4.98583, -4.89107, 1.2734],
                          [5.12371, -6.58952, 2.51502], [4.44504, -6.27396, 3.11854], [6.32436, -8.95849, 1.54133], [4.94189, -8.09755, 2.35112],
                          [6.48806, -6.24572, 3.11851], [7.26641, -5.53031, 2.44434], [6.73817, -6.71910, 4.25407], [4.09052, -8.26497, 1.80625],
                          [4.80079, -8.51063, 3.27822], [6.23058, -8.37401, 0.33742], [7.12646, -2.84192, 1.85399], [6.44668, -2.65528, 1.19010],
                          [7.84635, -2.49324, 1.60226], [7.32039, -3.79696, 1.89821], [6.73612, -2.33452, 3.19812], [7.44308, -2.46895, 3.79819],
                          [6.22916, -0.01794, 4.68140], [6.44309, -0.83751, 3.08529], [5.50488, -3.11338, 3.68042], [4.67676, -3.43939, 2.78112],
                          [5.39170, -3.34123, 4.89723], [5.61047, -0.71047, 2.54217], [7.19689, -0.39642, 2.59351], [7.57250, 0.05966, 4.83832]])
        chain = chain/10

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


