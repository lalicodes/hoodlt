"""
:module: sc_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the sc lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2015
"""


import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatSC(LatSimple):
    """Implementation of the simple cubic lattice

    Space Group: Pm-3m (221)

    Strukturbericht: Ah

    Pearson symbol: cP1
    """

    def __init__(self, l_value, a_nn_e):

        """The constructor

        :param l_value: linear size
        :param a_nn_e:
        """

        super(LatSC, self).__init__(l_value, a_nn_e)

        self.a_vector[0] = self.a_nn*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.a_nn*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.a_nn*np.array([0.0, 0.0, 1.0])

        # number of primitive vectors
        self.typ = np.array([1])
        # basis vectors
        self.v_vector = np.zeros([np.sum(self.typ), 3])
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])

    @staticmethod
    def name():
        """Lattice name

        :return: [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['sc', 'sc', '221', 'Ah']
