"""
:module: LCysteine
:platform: Unix, Windows
:synopsis: Implements a class defining Ethane

.. moduleauthor:: Jacob Austin <wjaustin2@iastate.edu> October 2019
.. history:
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs

class LCysteine(SolventAbs):
    """
    Defines solvent LCysteine  (:math:`\\mbox{CH}_3\\mbox{CH}_3`)
    """

    def __init__(self, ff):
        super(LCysteine, self).__init__(1, ff, 13, 'lcysteine')

        # atom types
        #self.types = ['C_l0', 'C_l1', 'C_l2', 'O_l1', 'O_l2', 'N_l', 'S_l',
        #              'H_l1', 'H_l3', 'H_lA', 'H_lB']
        #self.typeid = ['S_l','C_l0','C_l1','H_l1','C_l2','O_l1','O_l2','N_l','H_l3','H_l3','H_lA',
        #               'H_lB','H_lB']


        self.types = ['C_l0', 'C_l1', 'C_l2', 'O_l1', 'O_l2', 'N_l', 'S',
                      'H_l1', 'H_l3', 'H_lA', 'H_lB']
        self.typeid = ['S','C_l0','C_l1','H_l1','C_l2','O_l1','O_l2','N_l','H_l3','H_l3','H_lA',
                       'H_lB','H_lB']

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('S_l')
        self.mass[1] = self.ff_reader.get_molecular_weight('C_l0')
        self.mass[2] = self.ff_reader.get_molecular_weight('C_l1')
        self.mass[3] = self.ff_reader.get_molecular_weight('H_l1')
        self.mass[4] = self.ff_reader.get_molecular_weight('C_l2')
        self.mass[5] = self.ff_reader.get_molecular_weight('O_l1')
        self.mass[6] = self.ff_reader.get_molecular_weight('O_l2')
        self.mass[7] = self.ff_reader.get_molecular_weight('N_l')
        self.mass[8] = self.ff_reader.get_molecular_weight('H_l3')
        self.mass[9] = self.ff_reader.get_molecular_weight('H_l3')
        self.mass[10] = self.ff_reader.get_molecular_weight('H_lA')
        self.mass[11] = self.ff_reader.get_molecular_weight('H_lB')
        self.mass[12] = self.ff_reader.get_molecular_weight('H_lB')

        # charges
        self.charge[0] = self.ff_reader.get_charge('S_l')
        self.charge[1] = self.ff_reader.get_charge('C_l0')
        self.charge[2] = self.ff_reader.get_charge('C_l1')
        self.charge[3] = self.ff_reader.get_charge('H_l1')
        self.charge[4] = self.ff_reader.get_charge('C_l2')
        self.charge[5] = self.ff_reader.get_charge('O_l1')
        self.charge[6] = self.ff_reader.get_charge('O_l2')
        self.charge[7] = self.ff_reader.get_charge('N_l')
        self.charge[8] = self.ff_reader.get_charge('H_l3')
        self.charge[9] = self.ff_reader.get_charge('H_l3')
        self.charge[10] = self.ff_reader.get_charge('H_lA')
        self.charge[11] = self.ff_reader.get_charge('H_lB')
        self.charge[12] = self.ff_reader.get_charge('H_lB')


        # atom positions (nm)
        self.position = self.build_chain()

    def build_chain(self):
        """
        called by initializer to get the positions

        :return: the position of lcysteine
        """
        chain = np.array([[0,0,0], [-.0672,.1466,-.085],[-.2206,.1466,-.085],
                          [-.2591,.1466,.0175],[-.2733,.2745,-.152],
                          [-.22,.384,-.147],[-.3922, .2515,-.213],[-.2802,.0291,-.157],
                          [-.0283,.2359,-.035],[-.0288,.1472,-.188],[-.4115,.1551,-.202],
                          [-.244,-.055,-.11],[-.237,.0235,-.25]])
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
        bond = np.array([[0,1], [1,8],[1,9],[1,2],[2,3],[2,7],[7,11],[7,12],[2,4]
                         ,[4,5],[4,6],[6,10]])
        return bond

    def get_bond_types(self):
        """
        :return:
        """
        return ['S_l-C_l0','C_l0-H_l3','C_l0-H_l3','C_l0-C_l1','C_l1-H_l1','C_l1-N_l','N_l-H_lB',
                'N_l-H_lB','C_l1-C_l2','C_l2-O_l1','C_l2-O_l2','O_l2-H_lA']

    def get_angles(self):
        """ overrides method in LigandAbs
        :return: angles
        """
        angles = np.array([[0,1,8],[0,1,9],[0,1,2],[8,1,9],[8,1,2],[9,1,2],[1,2,3],[1,2,4],[3,2,7],
                            [3,2,4],[11,7,2],[12,7,2],[11,7,12],[7,2,4],[2,4,5],[2,4,6],[5,4,6],[4,6,10],[1,2,7]])
        return angles

    def get_angle_types(self):
        """
        :return:
        """
        return ['S_l-C_l0-H_l3','S_l-C_l0-H_l3','S_l-C_l0-C_l1','H_l3-C_l0-H_l3','H_l3-C_l0-C_l1','H_l3-C_l0-C_l1',
                'C_l0-C_l1-H_l1','C_l0-C_l1-C_l2','H_l1-C_l1-N_l','H_l1-C_l1-C_l2','H_lB-N_l-C_l1','H_lB-N_l-C_l1',
                'H_lB-N_l-H_lB','N_l-C_l1-C_l2','C_l1-C_l2-O_l1','C_l1-C_l2-O_l2','O_l1-C_l2-O_l2','C_l2-O_l2-H_lA','C_l0-C_l1-N_l']

    def get_dihedrals(self):
        """ overrides method in LigandAbs
        :return: dihedrals
        """
        dihedrals = np.array([[11,7,2,4],[12,7,2,4],[12,7,2,1],[11,7,2,1],[12,7,2,3],[11,7,2,3],[8,1,2,4],[9,1,2,4],
                              [8,1,2,3],[9,1,2,3],[8,1,2,7],[9,1,2,7],[10,6,4,2],[10,6,4,5],[7,2,4,6],[6,4,2,1],[0,1,2,4],[0,1,2,3],[0,1,2,7]])
        return dihedrals

    def get_dihedral_types(self):
        """

        :return:
        """
        return ['H_lB-N_l-C_l1-C_l2','H_lB-N_l-C_l1-C_l2','H_lB-N_l-C_l1-C_l0','H_lB-N_l-C_l1-C_l0','H_lB-N_l-C_l1-H_l1','H_lB-N_l-C_l1-H_l1',
                'H_l3-C_l0-C_l1-C_l2','H_l3-C_l0-C_l1-C_l2','H_l3-C_l0-C_l1-H_l1','H_l3-C_l0-C_l1-H_l1','H_l3-C_l0-C_l1-N_l','H_l3-C_l0-C_l1-N_l',
                'H_lA-O_l2-C_l2-C_l1','H_lA-O_l2-C_l2-O_l1','N_l-C_l1-C_l2-O_l2',
                'O_l2-C_l2-C_l1-C_l0','S_l-C_l0-C_l1-C_l2','S_l-C_l0-C_l1-H_l1','S_l-C_l0-C_l1-N_l']

    def get_impropers(self):
        """ overrides method in LigandAbs
        :return:
        """
        impropers = np.array([[6,4,2,5]])
        return impropers

    def get_improper_types(self):

        return ['O_l2-C_l2-C_l1-O_l1']

