"""
:module: chain_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the class implementing the chain lattice

.. moduleauthor:: Xun Zha <xzha@iastate.edu> November 2018
"""

import numpy as np
from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatChain(LatSimple):
    """Implementation of the chain lattice

    Space Group: -

    Strukturbericht: -

    Pearson symbol: -
    """

    def __init__(self, l_value, a_nn_e):
        """The constructor

        :param l_value: linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatChain, self).__init__(l_value, a_nn_e)

        self.l[1] = 1
        self.l[2] = 1

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis*np.array([0.0, 0.0, 1.0])

        # number of primitive vectors of type A
        self.typ = np.array([1])
        # basis vectors
        self.v_vector = np.zeros([np.sum(self.typ), 3])
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])

    def resize_box(self, axis_l):
        """resize the box in y and z dimension

        :param axis_l: new y- and z-axis length
        :return:
        """
        self.a_vector[1] *= axis_l/self.b_axis
        self.a_vector[2] *= axis_l/self.c_axis
        self.b_axis = axis_l
        self.c_axis = axis_l

    @staticmethod
    def name():
        """Lattice name

        :return: strings as [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['chain', 'chain', '-', '-']