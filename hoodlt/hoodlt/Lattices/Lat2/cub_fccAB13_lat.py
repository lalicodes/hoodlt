"""
:module: cub_fccAB13_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{cub-fccAB}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from numpy import array
from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary

class LatCubfccAB13Base56(LatBinary):
    """Implementation of the :math:`\\mbox{cub-fccAB}_{13}` lattice

    Fifty six particles per unit cell
    
    Space Group: Fm-3m (225)

    Wyckoff positions: A- 4a B- 48i  4b

    Strukturbericht: :math:`\\mbox{D}2_f`

    Pearson symbol: cF56"""

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCubfccAB13Base56, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, 1-np.sqrt(2.0/3.0), (np.sqrt(10)+1)/9.0, 1.0])

        if gamma < 1-np.sqrt(2.0/3.0):
            c_fac = np.sqrt(2)
            y_c = 0.5-1.0/np.sqrt(6.0)
        elif gamma > (np.sqrt(10)+1)/9.0:
            c_fac = 6*gamma/np.sqrt(2.0)
            y_c = 1/6.0
        else:
            y_c = (-2*gamma+np.sqrt(-2.0*gamma**2+4*gamma+2))*gamma/(4*gamma+2-6*gamma**2)
            c_fac = gamma/(np.sqrt(2)*y_c)

        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 4
        self.typ[1] = 52

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        norm = self.a_axis

        self.v_vector[0] = norm * array([0.0, 0.0, 0.0])
        self.v_vector[1] = norm * array([0.5, 0.5, 0.0])
        self.v_vector[2] = norm * array([0.0, 0.5, 0.5])
        self.v_vector[3] = norm * array([0.5, 0.0, 0.5])

        self.v_vector[4] = norm * array([0.5, y_c, y_c])
        self.v_vector[5] = norm * array([y_c, 0.5, y_c])
        self.v_vector[6] = norm * array([y_c, y_c, 0.5])

        self.v_vector[7] = norm * array([0.5, 1 - y_c, y_c])
        self.v_vector[8] = norm * array([y_c, 0.5, 1 - y_c])
        self.v_vector[9] = norm * array([1 - y_c, y_c, 0.5])

        self.v_vector[10] = norm * array([0.5, y_c, 1 - y_c])
        self.v_vector[11] = norm * array([1 - y_c, 0.5, y_c])
        self.v_vector[12] = norm * array([y_c, 1 - y_c, 0.5])

        self.v_vector[13] = norm * array([0.5, 1 - y_c, 1 - y_c])
        self.v_vector[14] = norm * array([1 - y_c, 0.5, 1 - y_c])
        self.v_vector[15] = norm * array([1 - y_c, 1 - y_c, 0.5])

        self.v_vector[16] = norm * (array([0.5, y_c, y_c]) + array([0.0, 0.5, 0.5]))
        self.v_vector[17] = norm * (array([y_c, 0.5, y_c]) + array([0.0, -0.5, 0.5]))
        self.v_vector[18] = norm * (array([y_c, y_c, 0.5]) + array([0.0, 0.5, -0.5]))
        self.v_vector[19] = norm * (array([0.5, 1 - y_c, y_c]) + array([0.0, -0.5, 0.5]))
        self.v_vector[20] = norm * (array([y_c, 0.5, 1 - y_c]) + array([0.0, -0.5, -0.5]))
        self.v_vector[21] = norm * (array([1 - y_c, y_c, 0.5]) + array([0.0, 0.5, -0.5]))
        self.v_vector[22] = norm * (array([0.5, y_c, 1 - y_c]) + array([0.0, 0.5, -0.5]))
        self.v_vector[23] = norm * (array([1 - y_c, 0.5, y_c]) + array([0.0, -0.5, 0.5]))
        self.v_vector[24] = norm * (array([y_c, 1 - y_c, 0.5]) + array([0.0, -0.5, -0.5]))
        self.v_vector[25] = norm * (array([0.5, 1 - y_c, 1 - y_c]) + array([0.0, -0.5, -0.5]))
        self.v_vector[26] = norm * (array([1 - y_c, 0.5, 1 - y_c]) + array([0.0, -0.5, -0.5]))
        self.v_vector[27] = norm * (array([1 - y_c, 1 - y_c, 0.5]) + array([0.0, -0.5, -0.5]))

        self.v_vector[28] = norm * (array([0.5, y_c, y_c]) + array([-0.5, 0.0, 0.5]))
        self.v_vector[29] = norm * (array([y_c, 0.5, y_c]) + array([0.5, 0.0, 0.5]))
        self.v_vector[30] = norm * (array([y_c, y_c, 0.5]) + array([0.5, 0.0, -0.5]))
        self.v_vector[31] = norm * (array([0.5, 1 - y_c, y_c]) + array([-0.5, 0.0, 0.5]))
        self.v_vector[32] = norm * (array([y_c, 0.5, 1 - y_c]) + array([0.5, 0.0, -0.5]))
        self.v_vector[33] = norm * (array([1 - y_c, y_c, 0.5]) + array([-0.5, 0.0, -0.5]))
        self.v_vector[34] = norm * (array([0.5, y_c, 1 - y_c]) + array([-0.5, 0.0, -0.5]))
        self.v_vector[35] = norm * (array([1 - y_c, 0.5, y_c]) + array([-0.5, 0.0, 0.5]))
        self.v_vector[36] = norm * (array([y_c, 1 - y_c, 0.5]) + array([0.5, 0.0, -0.5]))
        self.v_vector[37] = norm * (array([0.5, 1 - y_c, 1 - y_c]) + array([-0.5, 0.0, -0.5]))
        self.v_vector[38] = norm * (array([1 - y_c, 0.5, 1 - y_c]) + array([-0.5, 0.0, -0.5]))
        self.v_vector[39] = norm * (array([1 - y_c, 1 - y_c, 0.5]) + array([-0.5, 0.0, -0.5]))

        self.v_vector[40] = norm * (array([0.5, y_c, y_c]) + array([-0.5, 0.5, 0.0]))
        self.v_vector[41] = norm * (array([y_c, 0.5, y_c]) + array([0.5, -0.5, 0.0]))
        self.v_vector[42] = norm * (array([y_c, y_c, 0.5]) + array([0.5, 0.5, 0.0]))
        self.v_vector[43] = norm * (array([0.5, 1 - y_c, y_c]) + array([-0.5, -0.5, 0.0]))
        self.v_vector[44] = norm * (array([y_c, 0.5, 1 - y_c]) + array([0.5, -0.5, 0.0]))
        self.v_vector[45] = norm * (array([1 - y_c, y_c, 0.5]) + array([-0.5, 0.5, 0.0]))
        self.v_vector[46] = norm * (array([0.5, y_c, 1 - y_c]) + array([-0.5, 0.5, 0.0]))
        self.v_vector[47] = norm * (array([1 - y_c, 0.5, y_c]) + array([-0.5, -0.5, 0.0]))
        self.v_vector[48] = norm * (array([y_c, 1 - y_c, 0.5]) + array([0.5, -0.5, 0.0]))
        self.v_vector[49] = norm * (array([0.5, 1 - y_c, 1 - y_c]) + array([-0.5, -0.5, 0.0]))
        self.v_vector[50] = norm * (array([1 - y_c, 0.5, 1 - y_c]) + array([-0.5, -0.5, 0.0]))
        self.v_vector[51] = norm * (array([1 - y_c, 1 - y_c, 0.5]) + array([-0.5, -0.5, 0.0]))

        self.v_vector[52] = norm * array([0.5, 0.5, 0.5])
        self.v_vector[53] = norm * array([0.5, 0.0, 0.0])
        self.v_vector[54] = norm * array([0.0, 0.5, 0.0])
        self.v_vector[55] = norm * array([0.0, 0.0, 0.5])

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['cubfccAB13', 'cub-fccAB$_{13}$', '225', 'D2_f', '#FFFF00']
