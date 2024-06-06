"""
:module: LatticeSimple
:platform: Unix, Windows
:synopsis: Defines the abstract classes implementing binary lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
"""

from __future__ import division

import numpy as np

from hoodlt.D_matrix.Dmatrix_Mixt_Lattices import DMixtLattice


class LatSimple(DMixtLattice):
    """Implements a general simple lattice


    """

    def __init__(self, l_value, a_nn_e):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        """

        super(LatSimple, self).__init__(l_value, a_nn_e)
        self.radius[0] = 0.5*self.a_nn

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """

        return ['Simple', 'Simple', '-', '-']

    def pf(self):
        """Packing fraction

        :return: packing fraction
        :rtype: float
        """

        return np.pi*self.g_l()/6.0
