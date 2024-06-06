"""
:module: MgZn2_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{MgZn}_{2}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatMgZn2Base12(LatBinary):
    """Implementation of the :math:`\mbox{MgZn}_2` lattice

    Twelve particles per unit cell

    Space Group: :math:`P6_3/mmc` (194)

    Wyckoff positions: A- 4f B- 2a 6h

    Strukturbericht: C14

    Pearson symbol: hP12

    """
   
    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatMgZn2Base12, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, np.sqrt(2.0 / 3.0), 1.0])

        if gamma < np.sqrt(2.0/3.0):
            a_fac = 2.0*np.sqrt(2.0/3.0)
        else:
            a_fac = 2.0*gamma

        c_fac = np.sqrt(8.0/3.0)

        self.a_axis *= a_fac
        self.b_axis *= a_fac
        self.c_axis *= c_fac*a_fac

        # redefine self.l for hexagonal lattice, to make it more isotropic
        if isinstance(l_value, int):
            self.l[1] = int(np.rint(1.15*l_value))
            self.l[2] = int(np.rint(float(l_value*self.a_axis/self.c_axis)))

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([-0.5, 0.5*np.sqrt(3.0), 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 4
        self.typ[1] = 8

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = self.a_axis*np.array([0.0, 1.0/np.sqrt(3.0), c_fac/16.0])
        self.v_vector[1] = self.a_axis*np.array([0.0, 1.0/np.sqrt(3.0), c_fac*7.0/16.0])
        self.v_vector[2] = self.a_axis*np.array([1.0/2.0, 1.0/(2.0*np.sqrt(3.0)), c_fac*9.0/16.0])
        self.v_vector[3] = self.a_axis*np.array([1.0/2.0, 1.0/(2.0*np.sqrt(3.0)), c_fac*15.0/16.0])

        self.v_vector[4] = np.array([0.0, 0.0, 0.0])
        self.v_vector[5] = self.a_axis*np.array([1.0/2.0, 1.0/np.sqrt(3.0), c_fac/4.0])
        self.v_vector[6] = self.a_axis*np.array([1.0/4.0, np.sqrt(3.0)/12.0, c_fac/4.0])
        self.v_vector[7] = self.a_axis*np.array([3.0/4.0, np.sqrt(3.0)/12.0, c_fac/4.0])
        self.v_vector[8] = self.a_axis*np.array([0.0, 0.0, c_fac/2.0])
        self.v_vector[9] = self.a_axis*np.array([0.0, np.sqrt(3.0)/6.0, 3.0*c_fac/4.0])
        self.v_vector[10] = self.a_axis*np.array([-1.0/4.0, 5.0*np.sqrt(3.0)/12.0, 3.0*c_fac/4.0])
        self.v_vector[11] = self.a_axis*np.array([1.0/4.0, 5.0*np.sqrt(3.0)/12.0, 3.0*c_fac/4.0])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['MgZn2', 'MgZn$_2$', '194', 'C14', '#56A0D3']
