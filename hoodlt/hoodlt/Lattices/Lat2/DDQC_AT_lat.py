"""
:module: DDQC_AT_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{DDQC-AT}` 

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
"""

import numpy as np

from hoodlt.Lattices.Lat2.LatticeBinary import LatBinary


class LatDDQCAT(LatBinary):

    """Implementation of the :math:`\\mbox{DDQC-AT}` quasicrystal

    unit cell

    Space Group: --

    Wyckoff positions: --

    Strukturbericht: --

    Pearson symbol: --
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

        :param l_value: lattice size
        :param a_nn_e: spacing
        :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
        """

        super(LatDDQCAT, self).__init__(l_value, a_nn_e, gamma)

        self.gamma_crit = np.array([0.0, 0.414, 0.462, 1.0])

    def pf(self):
        """Redefine packing fraction

        :return: packing fraction
        :rtype: float
        """

        if self.gam <= self.gamma_crit[1]:
            val = 0.5*np.pi*(self.gam**3+1/6.0) + np.pi*(self.gam**3+0.5)/(3.0*np.sqrt(3))
        elif (self.gam < self.gamma_crit[2]) and (self.gam > self.gamma_crit[1]):
            val = (0.5*np.pi*(self.gam**3+1/6.0) + np.pi*(self.gam**3+0.5)/(3.0*np.sqrt(3)))/(self.gam*(1+np.sqrt(2)))
        else:
            val = 0.5*0.0463*(1.0/self.gam**3 + 6.0) + 0.5*0.0535*(1.0/self.gam**3 + 2.0)

        return val

    @staticmethod
    def name():
        """name

        :return: list of names
        :rtype: list
        """
        return ['DDQC', 'DDQC/AT', 'QC', 'QC', '#ED872D']
