"""
:module: Th3P4_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{Th}_3\\mbox{P}_4` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2019
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatTh3P4Base28(LatBinary):
    """Implements the lattice :math:`\\mbox{Th}_3\\mbox{P}_4` lattice

    twenty eight particles per unit cell

    Space Group: I43d (220)

    Wyckoff positions: A- 12b B- 16c

    Strukturbericht: D7_3

    Pearson symbol: cI28
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
        """

        super(LatTh3P4Base28, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(966)/21.0-1, 1.0])

        x_pos = 11.0/12.0

        if gamma < self.gamma_crit[1]:
            self.a_axis = 8*self.a_nn/np.sqrt(14)
        else:
            self.a_axis = 8*(self.a_nn/np.sqrt(14))*np.sqrt(21.0/46.0)*(1+gamma)

        # primitive vectors
        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.a_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.a_axis*np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 12
        self.typ[1] = 16

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))
        self.v_vector[0] = self.a_axis*np.array([7.0/8.0, 0.0, 0.25])
        self.v_vector[1] = self.a_axis*np.array([5.0/8.0, 0.0, 0.75])
        self.v_vector[2] = self.a_axis*np.array([0.25, 7.0/8.0, 0.0])
        self.v_vector[3] = self.a_axis*np.array([0.75, 5.0/8.0, 0.0])
        self.v_vector[4] = self.a_axis*np.array([0.0, 0.25, 7.0/8.0])
        self.v_vector[5] = self.a_axis*np.array([0.0, 0.75, 5.0/8.0])
        self.v_vector[6] = self.a_axis*np.array([3.0/8.0, 0.5, 0.75])
        self.v_vector[7] = self.a_axis*np.array([0.125, 0.5, 0.25])
        self.v_vector[8] = self.a_axis*np.array([0.75, 3.0/8.0, 0.5])
        self.v_vector[9] = self.a_axis*np.array([0.25, 0.125, 0.5])
        self.v_vector[10] = self.a_axis*np.array([0.5, 0.75, 3.0/8.0])
        self.v_vector[11] = self.a_axis*np.array([0.5, 0.25, 0.125])

        self.v_vector[12] = np.array([x_pos, x_pos, x_pos])
        self.v_vector[13] = np.array([-x_pos+0.5, -x_pos, x_pos+0.5])
        self.v_vector[14] = np.array([-x_pos, x_pos+0.5, -x_pos+0.5])
        self.v_vector[15] = np.array([x_pos+0.5, -x_pos+0.5, -x_pos])
        self.v_vector[16] = np.array([x_pos+0.25, x_pos+0.25, x_pos+0.25])
        self.v_vector[17] = np.array([-x_pos+0.25, -x_pos+0.75, x_pos+0.75])
        self.v_vector[18] = np.array([x_pos+0.75, -x_pos + 0.25, -x_pos+0.75])
        self.v_vector[19] = np.array([-x_pos+0.75, x_pos + 0.75, -x_pos + 0.25])

        self.v_vector[20:] = self.v_vector[12: 20]+ np.array([0.5, 0.5, 0.5])

        self.v_vector[12:] = self.a_axis*(self.v_vector[12:] % 1)

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """

        return ['Th3P4', 'Th$_3$P$_4$', '220', 'D7_3', '#D2691E', 'colour']
