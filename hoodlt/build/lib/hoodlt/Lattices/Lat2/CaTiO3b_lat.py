"""
:module: CaTiO3b_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{CaTiO}_3` lattices (binary version where Ti=O

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCaTiO3Base5(LatBinary):
    """Implementation of the :math:`\\mbox{CaTiO}_3` lattice
    
    Five particles per unit cell

    Space Group: Pm-3m (221)

    Wyckoff positions: A- 1a B- 1b 3c

    Strukturbericht: :math:`\mbox{E1}_2`

    Pearson symbol:  cP5
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCaTiO3Base5, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(2) - 1, 1/(2*np.sqrt(2)-1), 1.0])

        if gamma < self.gamma_crit[1]:
            c_fac = 1.0
        elif gamma < self.gamma_crit[2]:
            c_fac = (self.gam + 1.0) / np.sqrt(2)
        else:
            c_fac = 2*self.gam

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
        self.typ[1] = 4

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis * np.array([0.5, 0.5, 0.0])
        self.v_vector[2] = self.a_axis * np.array([0.5, 0.0, 0.5])
        self.v_vector[3] = self.a_axis * np.array([0.0, 0.5, 0.5])
        self.v_vector[4] = self.a_axis * np.array([0.5, 0.5, 0.5])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['CaTiO3 bin', 'CaTiO$_3$ (bin)', '221', 'E21', '#40E0D0']
