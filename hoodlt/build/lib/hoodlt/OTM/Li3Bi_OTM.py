"""
:module: MgZn2_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{Li}_{3}\\mbox{Bi}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.Li3Bi_lat import LatLi3BiBase16
import hoodlt.OTM.OTMobs as ObT


class OTMLatLi3BiBase16(LatLi3BiBase16):
    """ Implementation of the :math:`\\mbox{Li}_3\\mbox{Bi}` lattice
    sixteen particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, np.sqrt(2)-1.0, 1.0])
        super(OTMLatLi3BiBase16, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam < self.gamma_crit[1]:
            self.otm.indx = 0
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()

        elif self.gam < self.otm.gamma_crit[1]:
            eta_crit = np.pi*np.sqrt(3)*(1+3.0*self.gamma_crit[1]**3)/(4*(1+self.gamma_crit[1])**3)
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = eta_crit*(1.0+3.0*self.gam**3)/(1.0+3.0*self.gamma_crit[1]**3)
            self.otm.ratio_b = self.gamma_crit[1]/self.gam

        else:
            self.otm.indx = 1
            self.otm.gamma_bar = 0.5*np.sqrt(3)*(1+self.gam)-1
            self.otm.pf = self.pf()*8.0/(3.0*np.sqrt(3.0))
            self.otm.ratio_b = self.otm.gamma_bar/self.gam

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
