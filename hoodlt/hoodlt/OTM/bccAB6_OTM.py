"""
:module: bccAB6_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{bccAB}_{6}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.bccAB6_lat import LatBccAB6Base14
import hoodlt.OTM.OTMobs as ObT


class OTMLatbccAB6Base14(LatBccAB6Base14):
    """ Implementation of the :math:`\\mbox{bccAB}_6` lattice
    Fourteen particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, 1 / np.sqrt(6), np.sqrt(2)/(np.sqrt(10)-1), 1.0])
        super(OTMLatbccAB6Base14, self).__init__(l_value, a_nn_e, gamma)

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
            eta_crit = np.pi*np.sqrt(5)*5*(1+6.0*self.gamma_crit[1]**3)/(24*(1+self.gamma_crit[1])**3)
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = eta_crit*(1.0+6.0*self.gam**3)/(1.0+6.0*self.gamma_crit[1]**3)
            self.otm.ratio_b = self.gamma_crit[1]/self.gam

        elif self.gam < self.gamma_crit[2]:
            self.otm.indx = 1
            self.otm.gamma_bar = np.sqrt(10.0)*self.gam-1
            self.otm.pf = np.pi*np.sqrt(5)*5*(1+6*self.gam**3)/(24*(1+self.otm.gamma_bar)**3)
            self.otm.ratio_b = (np.sqrt(10)-1/self.gam)

        else:
            eta_crit = np.pi * np.sqrt(2) * (6 + 1.0 / self.gamma_crit[2] ** 3) / 96.0
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[2]
            self.otm.pf = eta_crit*(1.0+6.0*self.gam**3)/(1.0+6.0*self.gamma_crit[2]**3)
            self.otm.ratio_b = self.gamma_crit[2]/self.gam

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
