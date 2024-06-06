"""
:module: ReO3_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{ReO}_3` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatReO3Base4(LatBinary):
    """Implementation of the :math:`\\mbox{ReO}_3` lattice

    Four particles per unit cell
    
    Space Group: Pm-3m (221)

    Wyckoff positions: A- 1a B- 3c

    Strukturbericht: :math:`DO_9`

    Pearson symbol: cP4
    """
    
    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatReO3Base4, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, 1.0])

        c_fac = 1+gamma

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
        self.typ[1] = 3

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))
        
        self.v_vector[0] = self.a_axis*np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.5, 0.0, 0.0])
        self.v_vector[2] = self.a_axis*np.array([0.0, 0.5, 0.0])
        self.v_vector[3] = self.a_axis*np.array([0.0, 0.0, 0.5])
        
    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['ReO3', 'ReO$_3$', '224', 'DO9', '#808000']
