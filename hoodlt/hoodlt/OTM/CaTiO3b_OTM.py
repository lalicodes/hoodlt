"""
:module: CaTiO3b_OTM
:platform: Unix, Windows
:synopsis: Defines the class implementing the OTM :math:`\\mbox{CaTiO}_3` lattices (binary case)

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, Februrary 2020
"""

import numpy as np
from hoodlt.Lattices.Lat2.CaTiO3b_lat import LatCaTiO3Base5
import hoodlt.OTM.OTMobs as ObT


class OTMLatCaTiO3Base5(LatCaTiO3Base5):
    """ Implementation of the :math:`\\mbox{CaTiO}_3` lattice binary case
    five particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, 1.0/2.0, 1/(2*np.sqrt(2)-1), 1.0])
        super(OTMLatCaTiO3Base5, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables

        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam < self.gamma_crit[1]:
            self.otm.indx = 0
            self.otm.pf = self.pf()
            self.otm.gamma_bar = self.gam

        elif self.gam < self.otm.gamma_crit[1]:
            self.otm.indx = 1
            self.otm.pf = np.pi*(1+4*self.gam**3)/6
            self.otm.gamma_bar = self.gamma_crit[1]

        elif self.gam < self.gamma_crit[2]:
            self.otm.indx = 1
            self.otm.pf = np.pi*(4+1/self.gam**3)/48.0
            self.otm.gamma_bar = 2*np.sqrt(2)*self.gam - 1

        elif self.gam < self.gamma_crit[3]:
            self.otm.indx = 1
            self.otm.pf = np.pi*np.sqrt(2)*(1+4*self.gam**3)/(3*(1+self.gam)**3)
            self.otm.gamma_bar = (1+(1-np.sqrt(2))*self.gam)/np.sqrt(2)

        return self.otm
