"""
:module: bcc_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the bcc lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2015
"""

import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatBcc(LatSimple):
    """ Implementation of the bcc lattice with 1 base

    Space Group: Im-3m (229)

    Strukturbericht: A2

    Pearson symbol: cI2
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatBcc, self).__init__(l_value, a_nn_e)

        cfac = 2.0/np.sqrt(3)

        self.l[2] *= 2

        self.a_axis *= cfac
        self.b_axis *= cfac
        self.c_axis *= cfac

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.a_nn*np.array([1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)])

        # Number of primitive vectors of type A
        self.typ[0] = 1

        # basis vectors
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])

    def l_box(self):
        """ Returns the box size in the natural unit of length

        :return: lattice box in units
        :rtype: numpy.array
        """
        return self.l[0]*np.array([self.a_axis, self.b_axis, self.c_axis])

    @staticmethod
    def name():
        """Lattice name

        :return: [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['bcc', 'bcc', '229', 'A2']


class LatBccBase2(LatSimple):
    """ Implementation of the bcc lattice with 2 basis

    Space Group: Im-3m (229)

    Strukturbericht: A2

    Pearson symbol: cI2
    """

    def __init__(self, l_value, a_nn_e):

        """The constructor

        :param l_value: linear size
        :param a_nn_e: minimum A-A separation
        :return:
        """

        super(LatBccBase2, self).__init__(l_value, a_nn_e)

        cfac = 2/np.sqrt(3)

        self.a_axis *= cfac
        self.b_axis *= cfac
        self.c_axis *= cfac

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis*np.array([0.0, 0.0, 1.0])

        # Number of primitive vectors of type A
        self.typ[0] = 2
      
        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = 0.5*np.array([self.a_axis, self.b_axis, self.c_axis])

    @staticmethod
    def name():
        """Lattice name

        :return: [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['bcc', 'bcc', '229', 'A2']
