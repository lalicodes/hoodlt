"""
:module: hcp_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the hcp lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, September2014
"""

import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatHcpBase2(LatSimple):
    """Implementation of the hcp lattice

    Space Group: :math:`P6_3/mmc`

    Strukturbericht: A3

    Pearson symbol: hP2
    """

    def __init__(self, l_value, a_nn_e):
        """The constructor

        :param l_value: linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatHcpBase2, self).__init__(l_value, a_nn_e)

        if isinstance(l_value, int):
            l_x = l_value
            l_y = int(np.rint(1.15*l_value))
            l_z = int(np.rint(float(l_value*self.a_axis/self.c_axis)))
            self.l = np.array([l_x, l_y, l_z])

        self.c_axis *= 2*np.sqrt(2.0/3.0)

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([-0.5, 0.5*np.sqrt(3), 0.0])
        self.a_vector[2] = self.c_axis*np.array([0, 0, 1.0])

        # number of primitive vectors of type A
        self.typ = np.array([2])
      
        # basis vectors
        self.v_vector = np.zeros([np.sum(self.typ), 3])
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_nn*np.array([0.5, 0.5/np.sqrt(3.0), np.sqrt(2.0/3.0)])

    @staticmethod
    def name():
        """Lattice name

        :return: strings as [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['hcp', 'hcp', '194', 'A3']
