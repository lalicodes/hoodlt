"""
:module: AuCu_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{AuCu}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.AuCu_lat import LatAuCuBase4

import hoodlt.OTM.OTMobs as ObT


class OTMLatAuCuBase4(LatAuCuBase4):
    """ Implementation of the :math:`\\mbox{AuCu}` lattice
    four particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, 1.0/np.sqrt(2), 1.0])
        super(OTMLatAuCuBase4, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam > self.gamma_crit[1]:
            self.otm.indx = 0
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()

        elif self.gam > self.otm.gamma_crit[1]:
            self.otm.indx = 1
            g_b = (self.gam**2-np.sqrt(2.0*self.gam**4-self.gam**2))/(1-self.gam**2)
            self.otm.gamma_bar = g_b
            self.otm.pf = np.pi*(1+1/self.gam**3)*self.otm.gamma_bar**2*self.gam/6.0
            self.otm.ratio_a = 1.0/np.sqrt(g_b**2+2*g_b-1)

        else:
            self.otm.indx = 1
            self.otm.gamma_bar = 1.0
            self.otm.pf = 1.2
            self.otm.ratio_a = 1.0/np.sqrt(2)

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
