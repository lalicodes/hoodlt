"""
:module: Fe4C_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{CFe}_{4}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatCFe4Base5(LatBinary):
    """Implementation of the :math:`\\mbox{Fe}_4\\mbox{C}` lattice

    Five particles per unit cell
    
    Space Group: P-43m (215)

    Wyckoff positions: A- 1a B- 4e

    Strukturbericht: -

    Pearson symbol: cP5
    """
    
    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        super(LatCFe4Base5, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, (np.sqrt(3.0)-1)/(np.sqrt(0.5*3)+1), 1.0])

        if gamma > (np.sqrt(3.0)-1)/(np.sqrt(0.5*3)+1):
            u_val = 0.5*(1/gamma+1)/(np.sqrt(3.0/2.0)+1+1.0/gamma)
            c_fac = gamma/(np.sqrt(2.0)*(1-2*u_val))
        else:
            gamma_c = (np.sqrt(3)-1)/(np.sqrt(3.0/2.0)+1)
            u_val = 0.5*(1/gamma_c+1)/(np.sqrt(3.0/2.0)+1+1.0/gamma_c)
            c_fac = 1.0
            
        self.a_axis *= c_fac
        self.b_axis *= c_fac
        self.c_axis *= c_fac

        # primitive vectors
        self.a_vector[0] = self.a_axis * np.array([1.0, 0.0, 0.0])
        self.a_vector[1] = self.b_axis * np.array([0.0, 1.0, 0.0])
        self.a_vector[2] = self.c_axis * np.array([0.0, 0.0, 1.0])

        self.typ.resize(2)
        # Number of primitive vectors of type A and B
        self.typ[0] = 1
        self.typ[1] = 4

        # basis vectors
        self.v_vector = np.zeros((np.sum(self.typ), 3))

        self.v_vector[0] = np.array([0.0, 0.0, 0.0])
        self.v_vector[1] = self.a_axis*np.array([u_val, u_val, u_val])
        self.v_vector[2] = self.a_axis*np.array([1.0-u_val, 1.0-u_val, u_val])
        self.v_vector[3] = self.a_axis*np.array([1.0-u_val, u_val, 1.0-u_val])
        self.v_vector[4] = self.a_axis*np.array([u_val, 1.0-u_val, 1.0-u_val])
    
    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['Fe4C', 'Fe$_4$C', '215', '-', '#A7FC00']
