"""
:module: NaZn13_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{NaZn}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatNaZn13Base112(LatBinary):
    """Implementation of the :math:`\\mbox{NaZn}_{13}` lattice

    112 particles per unit cell

    Space Group: Fm-3c (226)

    Wyckoff positions: A- 8a B- 96i 8b

    Strukturbericht: :math:`\mbox{D}2_3`

    Pearson symbol: cF112
    """
    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatNaZn13Base112, self).__init__(l_value, a_nn_e, gamma)

        tau = 0.5 * (1 + np.sqrt(5.0))
        theta = 1 + 2 * (1 + tau) / (np.sqrt(1 + tau ** 2))
        psi = 2 * (11.0 + np.sqrt(5) + 2 * np.sqrt(10.0 * (1 + np.sqrt(5)))) / (5 + np.sqrt(5))
        gamma_c1 = (theta - np.sqrt(theta ** 2 - 6.0)) / 3.0
        gamma_c2 = (1 + np.sqrt(psi + 1)) / psi

        self.gamma_crit = np.array([0.0, gamma_c1, gamma_c2, 1.0])

        if self.gam < gamma_c1:
            c_fac = 2.0

        elif (self.gam >= gamma_c1) and (self.gam < gamma_c2):
            t_fac = (3.0 + np.sqrt(5))*np.sqrt(2.0/(5.0+np.sqrt(5.0)))
            c_fac = 2.0*(t_fac*self.gam + np.sqrt(3+6*self.gam+(8.0/np.sqrt(5)-5)*self.gam**2))/3.0

        else:
            c_fac = 2 * self.gam * (2 * tau + np.sqrt(tau ** 2 - 1)) / np.sqrt(1 + tau ** 2)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 8
        self.typ[1] = 104

        hb = np.zeros((8, 3))
        hb[1] = np.array([0.0, 0.5, 0.5])
        hb[2] = np.array([0.5, 0.5, 0.0])
        hb[3] = np.array([0.5, 0.0, 0.5])
        hb[4] = np.array([0.5, 0.5, 0.5])
        hb[5] = np.array([0.5, 0.0, 0.0])
        hb[6] = np.array([0.0, 0.5, 0.0])
        hb[7] = np.array([0.0, 0.0, 0.5])

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        # a vectors
        v_vec = np.array([0.25, 0.25, 0.25])
        for ind in range(8):
            self.v_vector[ind] = self.a_axis * (v_vec + hb[ind])

        modt = np.sqrt(1 + tau ** 2)
        # first icosahedra
        ico = np.zeros((13, 3))
        ico[1] = np.array([0, 1, tau])
        ico[2] = np.array([0, 1, -tau])
        ico[3] = np.array([0, -1, tau])
        ico[4] = np.array([0, -1, -tau])

        for ind in range(4):
            ico[5 + ind] = np.roll(ico[ind + 1], 1)
            ico[9 + ind] = np.roll(ico[ind + 1], 2)

        # translations of first icosahedra
        n_off = 8
        for ind1 in range(4):
            for ind2 in range(13):
                self.v_vector[n_off + ind2 + 13 * ind1] = a_nn_e * self.gam * ico[ind2] / modt + self.a_axis * hb[ind1]

        # second icosahedra
        ico = np.zeros((13, 3))
        ico[1] = np.array([1, 0, tau])
        ico[2] = np.array([-1, 0, tau])
        ico[3] = np.array([1, 0, -tau])
        ico[4] = np.array([-1, 0, -tau])

        for ind in range(4):
            ico[5 + ind] = np.roll(ico[ind + 1], 1)
            ico[9 + ind] = np.roll(ico[ind + 1], 2)

        # translations of second icosahedra
        n_off = 60
        for ind1 in range(4):
            for ind2 in range(13):
                self.v_vector[n_off + ind2 + 13 * ind1] = a_nn_e * self.gam * ico[ind2] / modt + self.a_axis * hb[4 + ind1]

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['NaZn13', 'NaZn$_{13}$', '226', 'D2_3', '#483C32']
