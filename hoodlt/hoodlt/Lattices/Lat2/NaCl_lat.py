"""
:module: NaCl_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{NaCl}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatNaClBase8(LatBinary):
    """Implementation of the :math:`\mbox{NaCl}` lattice

    Twelve particles per unit cell
    
    Space Group: Fm-3m (225)

    Wyckoff positions: A- 4a B- 4b

    Strukturbericht: B1

    Pearson symbol: cF8
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatNaClBase8, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(2.0)-1, 1.0])

        if gamma < np.sqrt(2.0)-1:
            c_fac = np.sqrt(2)
        else:
            c_fac = (1+gamma)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 4
        self.typ[1] = 4

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.5, 0.5, 0.0])
        self.v_vector[2] = self.a_axis*np.array([0.5, 0.0, 0.5])
        self.v_vector[3] = self.a_axis*np.array([0.0, 0.5, 0.5])
        self.v_vector[4] = self.a_axis*np.array([0.5, 0.0, 0.0])
        self.v_vector[5] = self.a_axis*np.array([0.0, 0.5, 0.0])
        self.v_vector[6] = self.a_axis*np.array([0.0, 0.0, 0.5])
        self.v_vector[7] = self.a_axis*np.array([0.5, 0.5, 0.5])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['NaCl', 'NaCl', '225', 'B1', '#FF2400']
