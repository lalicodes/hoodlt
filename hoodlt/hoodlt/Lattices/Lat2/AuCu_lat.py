"""
:module: AuCu_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{AuCu}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatAuCuBase4(LatBinary):

    """Implementation of the :math:`\\mbox{AuCu}` lattice

    Four particles per unit cell

    Space Group: P4/mmm (123)

    Wyckoff positions: A- 1a 1b B- 2e

    Strukturbericht: :math:`\mbox{L1}_0`

    Pearson symbol: cP4
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
        """

        super(LatAuCuBase4, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(3.0)-1, 1.0])
        
        if (0 < gamma) and (gamma < np.sqrt(3.0)-1):
            self.a_axis *= np.sqrt(2)
            self.b_axis *= np.sqrt(2)
        
        elif (gamma >= np.sqrt(3.0)-1) and (gamma <= 1.0):
            self.a_axis *= np.sqrt(2)
            self.b_axis *= np.sqrt(2)
            self.c_axis *= np.sqrt(2)*np.sqrt(0.5*gamma**2+gamma-0.5)

        else:
            raise ValueError('Parameter gamma is outside allowed range')

        # redefine self.l for hexagonal lattice, to make it more isotropic
        if isinstance(l_value, int):
            self.l[2] = int(np.rint(float(l_value*self.a_axis/self.c_axis)))

        # primitive vectors
        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis*np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 2
        self.typ[1] = 2

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.5, 0.5, 0.0])
        self.v_vector[2] = np.array([0.5*self.a_axis, 0.0, 0.5*self.c_axis])
        self.v_vector[3] = np.array([0.0, 0.5*self.a_axis, 0.5*self.c_axis])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['AuCu', 'AuCu', '123', 'L1_0', '#E30B5C']
