"""
:module: CU3Au_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{AuCu}_3` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.Cu3Au_lat import LatCu3AuBase4

import hoodlt.OTM.OTMobs as ObT


class OTMLatCu3AuBase4(LatCu3AuBase4):
    """ Implementation of the :math:`\\mbox{AuCu}_3` lattice
    three particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, 1/np.sqrt(2), 1.0])
        super(OTMLatCu3AuBase4, self).__init__(l_value, a_nn_e, gamma)

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
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = np.sqrt(2)*np.pi*(1+3*self.gam**3)/(3.0*(1+self.otm.gamma_bar)**3)
            self.otm.ratio_b = self.otm.gamma_bar/self.gam
        else:
            self.otm.indx = 1
            self.otm.gamma_bar = 2*self.gam-1
            self.otm.pf = np.sqrt(2)*np.pi*(1+3*self.gam**3)/(3.0*(1+self.otm.gamma_bar)**3)
            self.otm.ratio_b = self.otm.gamma_bar/self.gam

        if self.otm.pf > 1:
            self.otm.indx = -1
            
        return self.otm
