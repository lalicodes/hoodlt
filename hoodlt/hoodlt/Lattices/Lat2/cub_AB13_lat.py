"""
:module: cub_AB13_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{cubAB}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCubAB13Base14(LatBinary):
    """Implementation of the :math:`\\mbox{cub-AB}_{13}`

    Fourteen particles per unit cell
    
    Space Group: Pm-3m (221)

    Wyckoff positions: A- 1a B- 12i  1b

    Strukturbericht: -

    Pearson symbol: cP14
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCubAB13Base14, self).__init__(l_value, a_nn_e, gamma)

        g_crit_1 = (1 + 2.0 * np.sqrt(2)) / 3.0 - np.sqrt(3 + 4 * np.sqrt(2)) / 3.0
        g_crit_2 = 1.0 / (np.sqrt(5 + 2 * np.sqrt(2)) - 1)

        self.gamma_crit = np.array([0.0, g_crit_1, g_crit_2, 1.0])

        if gamma < g_crit_1:
            c_fac = 1.0
            x_f = g_crit_1 / np.sqrt(2)
        elif gamma > g_crit_2:
            c_fac = gamma * (1 + np.sqrt(2))
            x_f = 1 / (2 + np.sqrt(2))
        else:
            c_fac = np.sqrt(2) * (3 * gamma ** 2 - 2 * gamma - 1) / (
            4 * gamma - np.sqrt(-2 * gamma ** 2 + 12 * gamma + 6))
            x_f = 0.5 * gamma * (4 * gamma - np.sqrt(-2 * gamma ** 2 + 12 * gamma + 6)) / (
            3 * gamma ** 2 - 2 * gamma - 1)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 1
        self.typ[1] = 13

        norm = self.a_axis

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])

        self.v_vector[1] = norm * np.array([0.5, 0.5, 0.5])

        self.v_vector[2] = norm * (np.array([0.0, x_f, x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[3] = norm * (np.array([0.0, -x_f, x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[4] = norm * (np.array([0.0, -x_f, -x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[5] = norm * (np.array([0.0, x_f, -x_f]) + np.array([0.5, 0.5, 0.5]))

        self.v_vector[6] = norm * (np.array([x_f, 0.0, x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[7] = norm * (np.array([x_f, 0.0, -x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[8] = norm * (np.array([-x_f, 0.0, -x_f]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[9] = norm * (np.array([-x_f, 0.0, x_f]) + np.array([0.5, 0.5, 0.5]))

        self.v_vector[10] = norm * (np.array([x_f, x_f, 0.0]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[11] = norm * (np.array([-x_f, x_f, 0.0]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[12] = norm * (np.array([-x_f, -x_f, 0.0]) + np.array([0.5, 0.5, 0.5]))
        self.v_vector[13] = norm * (np.array([x_f, -x_f, 0.0]) + np.array([0.5, 0.5, 0.5]))

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['cubAB13', 'cub-AB$_{13}$', '221', '-', '#FF91A4']
