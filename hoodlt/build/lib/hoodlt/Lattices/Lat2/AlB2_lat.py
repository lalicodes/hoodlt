"""
:module: AlB2_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{AlB}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatAlB2Base3(LatBinary):
    """Implements the lattice :math:`\\mbox{AlB}_2`

    Three particles per unit cell

    Space Group: P6/mmm (191)

    Wyckoff positions: A- 1a B- 2d

    Strukturbericht: C32

    Pearson symbol: hP3
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
        """

        super(LatAlB2Base3, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(7.0/3.0)-1, 1/np.sqrt(3.0), np.sqrt(7.0/3.0)-1, 1.0])

        if gamma < np.sqrt(7.0/3.0)-1:
            self.a_axis = self.a_nn
            self.c_axis = self.a_nn
        elif gamma > 2.0/3.0:
            self.a_axis = np.sqrt(3.0)*gamma*self.a_nn
            self.c_axis = self.a_nn
        elif (gamma > 1/np.sqrt(3.0)) and (gamma <= 2.0/3.0):
            c_r = np.sqrt(-3.0*gamma**2+2.0*gamma+1.0)
            self.a_axis = np.sqrt(3.0)*gamma*self.a_nn
            self.c_axis = c_r*self.a_nn
        else:
            c_r = np.sqrt(gamma**2+2*gamma-1.0/3.0)
            self.a_axis = self.a_nn
            self.c_axis = c_r*self.a_nn

        self.b_axis = self.a_axis

        # redefine self.l for hexagonal lattice, to make it more isotropic
        if isinstance(l_value, int):
            self.l[1] = int(np.rint(1.15*l_value))
            self.l[2] = int(np.rint(float(l_value*self.a_axis/self.c_axis)))

        # primitive vectors
        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.a_axis*np.array([-0.5, 0.5*np.sqrt(3), 0.0])
        self.a_vector[2] = self.c_axis*np.array([0, 0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 1
        self.typ[1] = 2

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = np.array([0.5*self.a_axis, 0.5*self.a_axis/np.sqrt(3.0), 0.5*self.c_axis])
        self.v_vector[2] = np.array([0.0, self.a_axis/np.sqrt(3.0), 0.5*self.c_axis])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """

        return ['AlB2', 'AlB$_2$', '191', 'C32', '#002147', 'colour']
