"""
:module: A15_lat_Mixt
:platform: Unix, Windows
:synopsis: Defines the classes implementing the A15 lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2015
"""


import numpy as np

from hoodlt.Lattices.Lat1.LatticeSimple import LatSimple


class LatA15Base8(LatSimple):
    """ Implementation of the A15 lattice

    Space Group: Pm-3n (223)

    Strukturbericht: A15

    Pearson symbol: cP8
    """

    def __init__(self, l_value, a_nn_e):
        """ The constructor

        :param l_value: lattice linear size
        :param a_nn_e: minimum A-A separation
        """

        super(LatA15Base8, self).__init__(l_value, a_nn_e)

        c_fac = 2

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
        self.v_vector[1] = self.a_axis*np.array([0.5, 0.5, 0.5])
        self.v_vector[2] = self.a_axis*np.array([0.0, 0.25, 0.5])
        self.v_vector[3] = self.a_axis*np.array([0.0, 0.75, 0.5])
        self.v_vector[4] = self.a_axis*np.array([0.5, 0.0, 0.25])
        self.v_vector[5] = self.a_axis*np.array([0.5, 0.0, 0.75])
        self.v_vector[6] = self.a_axis*np.array([0.25, 0.5, 0.0])
        self.v_vector[7] = self.a_axis*np.array([0.75, 0.5, 0.0])


    @staticmethod
    def name():
        return ['A15', 'A15', '223', 'A15']
