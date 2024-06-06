"""
:module: bccAB6_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{CaB}_6` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCaB6Base7(LatBinary):
    """Implementation of the :math:`\\mbox{CaB}_6` lattice

    Seven particles per unit cell

    Space Group: Pm-3m (221)

    Wyckof Positions: A-1a B-6f

    Strukturbericht: :math:`D2_1`

    Pearson symbol: cP7

    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCaB6Base7, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, 1.0 / (1 + np.sqrt(2)), 1.0])

        if gamma < 1 / (1 + np.sqrt(2)):
            c_fac = 1
            u_val = (1 - np.sqrt(2) * gamma) / 2.0
        else:
            c_fac = gamma * (1 + np.sqrt(2))
            u_val = 1 / (2.0 * np.sqrt(2) + 2.0)

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
        self.typ[1] = 6

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis * np.array([0.5, 0.5, u_val])
        self.v_vector[2] = self.a_axis * np.array([0.5, 0.5, 1.0 - u_val])
        self.v_vector[3] = self.a_axis * np.array([u_val, 0.5, 0.5])
        self.v_vector[4] = self.a_axis * np.array([1 - u_val, 0.5, 0.5])
        self.v_vector[5] = self.a_axis * np.array([0.5, u_val, 0.5])
        self.v_vector[6] = self.a_axis * np.array([0.5, 1 - u_val, 0.5])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['CaB6', 'CaB$_6$', '221', 'D_21', '#3FFF00']
