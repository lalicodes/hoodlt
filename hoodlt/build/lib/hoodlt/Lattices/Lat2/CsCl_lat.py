"""
:module: CsCl_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{CsCl}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCsClBase2(LatBinary):

    """Implementation of the :math:`\mbox{CsCl}` lattice

    Two particles per unit cell

    Space Group: Pm-3m (221)

    Wyckoff positions: A- 1a B- 1b

    Strukturbericht: B2

    Pearson symbol:  cP2
    """
    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCsClBase2, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(3.0)-1, 1.0])

        if gamma < np.sqrt(3.0)-1:
            c_fac = 1.0
        else:
            c_fac = (1+gamma)/np.sqrt(3.0)

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
        self.typ[1] = 1

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = 0.5*self.a_axis*np.array([1.0, 1.0, 1.0])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['CsCl', 'CsCl', '221', 'B2', '#800000']
