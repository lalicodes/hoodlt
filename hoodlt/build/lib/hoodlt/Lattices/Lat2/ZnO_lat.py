"""
:module: ZnO_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{ZnO}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatZnOBase4(LatBinary):
    """Implementation of the :math:`\\mbox{ZnO}` lattice

    Eight particles per unit cell
    
    Space Group: P6_3mc (186)

    Wyckoff positions: A- 2a B- 2b

    Strukturbericht: B4

    Pearson symbol: hP4
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatZnOBase4, self).__init__(l_value, a_nn_e, gamma)

        # this condition imposes that both big and small particles form tetrahedra
        c_fac = np.sqrt(8.0/3.0)

        self.gam_crit = np.array([0.0, np.sqrt(3.0/2.0)-1, 1.0])

        if gamma > np.sqrt(3.0/2.0)-1:
            a_fac = 0.5*(1+gamma)*c_fac
            z_fac = 3.0/8.0
        else:
            a_fac = 1.0
            z_fac = 3.0/8.0

        # redefine self.l for hexagonal lattice, to make it more isotropic
        if isinstance(l_value, int):
            self.l[1] = int(np.rint(1.15*l_value))
            self.l[2] = int(np.rint(float(l_value*self.a_axis/self.c_axis)))

        self.a_axis *= a_fac
        self.b_axis *= a_fac
        self.c_axis *= c_fac*a_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis*np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.a_axis*np.array([-0.5, 0.5*np.sqrt(3), 0.0])
        self.a_vector[2] = self.c_axis*np.array([0, 0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 2
        self.typ[1] = 2

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([0.0, 1.0/np.sqrt(3.0), 0.5*c_fac])
        self.v_vector[2] = self.a_axis*np.array([0.0, 0.0, z_fac*c_fac])
        self.v_vector[3] = self.a_axis*np.array([0.0, 1.0/np.sqrt(3.0), (z_fac+0.5)*c_fac])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['ZnO', 'ZnO', '186', 'B4', '#C80815']
