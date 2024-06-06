"""
:module: fcc_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the fcc lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2015
"""

import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatFcc(LatSimple):
    """Implementation of the fcc lattice with one base

    Space Group: Fm-3m (225)

    Strukturbericht: A1

    Pearson symbol: cF4
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatFcc, self).__init__(l_value, a_nn_e)

        c_fac = np.sqrt(2)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        self.l[1] *= 2
        self.l[2] *= 2

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.a_nn*np.array([1/np.sqrt(2.0), 1/np.sqrt(2.0), 0.0])
        self.a_vector[2] = self.a_nn*np.array([1.0/np.sqrt(2.0), 0.0, 1.0/np.sqrt(2.0)])

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
        return ['fcc', 'fcc', '225', 'A1']


class LatFccBase4(LatSimple):
    """Implementation of the fcc lattice with four basis

    Space Group: Fm-3m (225)

    Strukturbericht: A1

    Pearson symbol: cF4
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatFccBase4, self).__init__(l_value, a_nn_e)

        cfac = np.sqrt(2)

        self.a_axis *= cfac
        self.b_axis *= cfac
        self.c_axis *= cfac

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis*np.array([0.0, 0.0, 1.0])

        self.typ[0] = 4

        self.v_vector = np.zeros((np.sum(self.typ), 3))
      
        # basis vectors
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_nn*np.array([0.0, 1/np.sqrt(2.0), 1/np.sqrt(2.0)])
        self.v_vector[2] = self.a_nn*np.array([1/np.sqrt(2.0), 0.0, 1/np.sqrt(2.0)])
        self.v_vector[3] = self.a_nn*np.array([1/np.sqrt(2.0), 1/np.sqrt(2.0), 0.0])

    @staticmethod
    def name():
        """Lattice name

        :return: [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['fcc', 'fcc', '225', 'A1']
