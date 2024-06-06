"""
:module: MgCu2_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{MgCu}_{2}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatMgCu2Base24(LatBinary):
    """ Implementation of the :math:`\\mbox{MgCu}_2` lattice

    Twenty four particles per unit cell
    
    Space Group: Fd-3m (227)

    Wyckoff positions: A- 8a B- 16d

    Strukturbericht: C14

    Pearson symbol: cF24
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatMgCu2Base24, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(2.0 / 3.0), 1.0])

        if gamma < np.sqrt(2.0 / 3.0):
            c_fac = 4.0 / np.sqrt(3.0)
        else:
            c_fac = 4.0 * gamma / np.sqrt(2.0)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 8
        self.typ[1] = 16

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = self.a_axis * np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis * np.array([0.0, 0.5, 0.5])
        self.v_vector[2] = self.a_axis * np.array([0.5, 0.0, 0.5])
        self.v_vector[3] = self.a_axis * np.array([0.5, 0.5, 0.0])
        self.v_vector[4] = self.a_axis * np.array([0.25, 0.25, 0.25])
        self.v_vector[5] = self.a_axis * np.array([0.25, 0.75, 0.75])
        self.v_vector[6] = self.a_axis * np.array([0.75, 0.25, 0.75])
        self.v_vector[7] = self.a_axis * np.array([0.75, 0.75, 0.25])

        self.v_vector[8] = self.a_axis * np.array([5.0/8.0, 5.0/8.0, 5.0/8.0])
        self.v_vector[9] = self.a_axis * np.array([5.0/8.0, 3.0/8.0, 3.0/8.0])
        self.v_vector[10] = self.a_axis * np.array([3.0/8.0, 5.0/8.0, 3.0/8.0])
        self.v_vector[11] = self.a_axis * np.array([3.0/8.0, 3.0/8.0, 5.0/8.0])
        self.v_vector[12] = self.a_axis * np.array([5.0/8.0, 1.0/8.0, 1.0/8.0])
        self.v_vector[13] = self.a_axis * np.array([5.0/8.0, 7.0/8.0, 7.0/8.0])
        self.v_vector[14] = self.a_axis * np.array([3.0/8.0, 1.0/8.0, 7.0/8.0])
        self.v_vector[15] = self.a_axis * np.array([3.0/8.0, 7.0/8.0, 1.0/8.0])
        self.v_vector[16] = self.a_axis * np.array([1.0/8.0, 5.0/8.0, 1.0/8.0])
        self.v_vector[17] = self.a_axis * np.array([1.0/8.0, 3.0/8.0, 7.0/8.0])
        self.v_vector[18] = self.a_axis * np.array([7.0/8.0, 5.0/8.0, 7.0/8.0])
        self.v_vector[19] = self.a_axis * np.array([7.0/8.0, 3.0/8.0, 1.0/8.0])
        self.v_vector[20] = self.a_axis * np.array([1.0/8.0, 1.0/8.0, 5.0/8.0])
        self.v_vector[21] = self.a_axis * np.array([1.0/8.0, 7.0/8.0, 3.0/8.0])
        self.v_vector[22] = self.a_axis * np.array([7.0/8.0, 1.0/8.0, 3.0/8.0])
        self.v_vector[23] = self.a_axis * np.array([7.0/8.0, 7.0/8.0, 5.0/8.0])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['MgCu2', 'MgCu$_2$', '227', 'C14', '#56A0D3']
