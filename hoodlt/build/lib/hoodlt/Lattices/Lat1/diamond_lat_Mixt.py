"""
:module: diamond_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the diamond lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2015
"""


import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatDiamondBase2(LatSimple):
    """Implementation of the diamond lattice with 2 basis

    Space Group: Fd-3m (227)

    Strukturbericht: A4

    Pearson symbol: cF8
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatDiamondBase2, self).__init__(l_value, a_nn_e)

        self.l[1] *= 2
        self.l[2] *= 2

        c_fac = 4.0/np.sqrt(3)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])

        self.a_vector[1] = self.c_axis*np.array([0.5, 0.5, 0.0])

        self.a_vector[2] = self.b_axis*np.array([0.5, 0.0, 0.5])
                       
        # Number of primitive vectors of type A
        self.typ[0] = 2

        self.v_vector = np.zeros((2, 3))
      
        # basis vectors
        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.25, 0.25, 0.25])

    def l_box(self):
        """ Returns the box size in the natural unit of length

        :return: lattice box in units
        :rtype: numpy.array
        """
        return self.l[0]*np.array([self.a_axis, self.b_axis, self.c_axis])

    @staticmethod
    def name():
        return ['Diamond', 'Diamond', '227', 'A4']


class LatDiamondBase8(LatSimple):
    """Implementation of the diamond lattice with 8 basis

    Space Group: Fd-3m (227)

    Strukturbericht: A4

    Pearson symbol: cF8
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatDiamondBase8, self).__init__(l_value, a_nn_e)

        c_fac = 4.0/np.sqrt(3)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis*np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis*np.array([0.0, 0.0, 1.0])

        # Number of primitive vectors of type A
        self.typ[0] = 8
      
        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.0, 0.5, 0.5])
        self.v_vector[2] = self.a_axis*np.array([0.5, 0.0, 0.5])
        self.v_vector[3] = self.a_axis*np.array([0.5, 0.5, 0.0])
        self.v_vector[4] = self.a_axis*np.array([0.25, 0.25, 0.25])
        self.v_vector[5] = self.a_axis*np.array([0.25, 0.75, 0.75])
        self.v_vector[6] = self.a_axis*np.array([0.75, 0.25, 0.75])
        self.v_vector[7] = self.a_axis*np.array([0.75, 0.75, 0.25])
     
    @staticmethod
    def name():
        """Lattice name

        :return: [name, latex name, space group, Strukturbericht]
        :rtype: list
        """
        return ['Diamond', 'Diamond', '227', 'A4']
