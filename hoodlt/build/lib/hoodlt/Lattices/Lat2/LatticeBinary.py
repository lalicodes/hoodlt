"""
:module: LatticeBinary
:platform: Unix, Windows
:synopsis: Defines the abstract classes implementing binary lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
"""

from __future__ import division

import numpy as np

from hoodlt.D_matrix.Dmatrix_Mixt_Lattices import DMixtLattice


class LatBinary(DMixtLattice):
    """Implements a general binary lattice


    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
        """

        super(LatBinary, self).__init__(l_value, a_nn_e)

        # math:`\\gamma`gamma value
        self.gam = gamma

        # abstract value for the critical :math:`\\gamma`
        self.gamma_crit = np.array([0.0, 1.0])

        self.radius.resize(2)
        self.radius[0] = 0.5*self.a_nn
        self.radius[1] = 0.5*self.a_nn*self.gam

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """

        return ['Binary', 'Binary', '-', '-', '#0000']

    def pf(self):
        """Packing fraction

        :return: packing fraction
        :rtype: float
        """

        return np.pi*self.g_l()*(self.typ[0] + self.typ[1]*self.gam**3)/(6.0*np.sum(self.typ))
