"""
:module: CaCu5_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{CaCu}_5` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCaCu5Base6(LatBinary):
    """Implementation of the :math:`\\mbox{CaCu}_5` lattice

    Six particles per unit cell

    Space Group: P6/mmm (191)

    Wyckof Positions: A-1a B-2c 3g

    Strukturbericht: :math:`\\mbox{D2}_d`

    Pearson symbol: hP6
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCaCu5Base6, self).__init__(l_value, a_nn_e, gamma)

        # gamma-crit
        self.gamma_crit = np.array([0.0, 2.0/np.sqrt(3.0)-1, (1+2*np.sqrt(19))/15.0, 1/(4.0/np.sqrt(3) - 1.0), 1.0])

        if gamma < 2.0 / np.sqrt(3.0) - 1:
            a_fac = 1.0
            c_fac = 1.0
        elif (gamma >= 2.0 / np.sqrt(3.0) - 1) and (gamma < (1 + 2 * np.sqrt(19)) / 15.0):
            a_fac = np.sqrt(3) * (1 + gamma) * 0.5
            c_fac = 1.0
        elif (gamma >= (1 + 2 * np.sqrt(19)) / 15.0) and (gamma < 1 / (4.0 / np.sqrt(3) - 1.0)):
            a_fac = np.sqrt(3) * (1 + gamma) * 0.5
            c_fac = 0.5 * np.sqrt(15.0 * gamma ** 2 - 2 * gamma - 1)
        else:
            a_fac = 2 * gamma
            c_fac = np.sqrt(8.0 / 3.0) * gamma

            # redefine self.l for hexagonal lattice, to make it more isotropic
        if isinstance(l_value, int):
            self.l[1] = int(np.rint(1.15 * l_value))
            self.l[2] = int(np.rint(float(l_value * self.a_axis / self.c_axis)))

        self.a_axis *= a_fac
        self.b_axis *= a_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([-0.5, 0.5 * np.sqrt(3), 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 1
        self.typ[1] = 5

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = np.array([0.5 * self.a_axis, 0.5 * self.a_axis / np.sqrt(3.0), 0.0])
        self.v_vector[2] = np.array([0.0, self.a_axis / np.sqrt(3.0), 0.0])
        self.v_vector[3] = np.array([0.5 * self.a_axis, 0.0, 0.5 * self.c_axis])
        self.v_vector[4] = np.array([0.25 * self.a_axis, np.sqrt(3.0) * self.a_axis / 4.0, 0.5 * self.c_axis])
        self.v_vector[5] = np.array([0.75 * self.a_axis, np.sqrt(3.0) * self.a_axis / 4.0, 0.5 * self.c_axis])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['CaCu5', 'CaCu$_5$', '191', 'D_2d', '#FF1DCE']
