"""
:module: Cr3Si_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{Cr}_3\\mbox{Si}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary

class LatSiCr3Base8(LatBinary):
    """Implementation of the :math:`\\mbox{Cr}_3\\mbox{Si}`

    Eight particles per unit cell

    Space Group: Pm-3n (223)

    Wyckof Positions: A-2a B-6c

    Strukturbericht: A15

    Pearson symbol: cP8
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatSiCr3Base8, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(5.0/3.0)-1, 1.0/(np.sqrt(5.0)-1.0), 1.0])

        u_val = 0.25

        if gamma > 1.0/(np.sqrt(5.0)-1.0):
            c_fac = 2*gamma
        elif gamma < np.sqrt(5.0/3.0)-1:
            c_fac = 2/np.sqrt(3.0)
        else:
            c_fac = 2*(1+gamma)/np.sqrt(5.0)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 2
        self.typ[1] = 6

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.5, 0.5, 0.5])
        self.v_vector[2] = self.a_axis*np.array([u_val, 0.0, 0.5])
        self.v_vector[3] = self.a_axis*np.array([1-u_val, 0.0, 0.5])
        self.v_vector[4] = self.a_axis*np.array([0.5, u_val, 0.0])
        self.v_vector[5] = self.a_axis*np.array([0.5, 1-u_val, 0.0])
        self.v_vector[6] = self.a_axis*np.array([0.0, 0.5, u_val])
        self.v_vector[7] = self.a_axis*np.array([0.0, 0.5, 1-u_val])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['Cr3Si', 'Cr$_3$Si', '223', 'A15', '#0000CD']
