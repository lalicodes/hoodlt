"""
:module: 
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{MgZn}_{2}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.MgZn2_lat import LatMgZn2Base12
import hoodlt.OTM.OTMobs as ObT


class OTMLatMgZn2Base12(LatMgZn2Base12):
    """ Implementation of the :math:`\\mbox{MgZn}_2` lattice
    Twelve particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, np.sqrt(2)/(np.sqrt(11)-np.sqrt(2)), 1.0])
        super(OTMLatMgZn2Base12, self).__init__(l_value, a_nn_e, gamma)

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

        elif self.gam >= self.otm.gamma_crit[1]:
            eta_crit = np.pi*np.sqrt(3)*(1+2.0*self.gamma_crit[1]**3)/16.0
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = eta_crit*(1.0/self.gam**3+2.0)/(1.0/self.gamma_crit[1]**3+2.0)
            self.otm.ratio_a = self.gam/self.gamma_crit[1]

        else:
            self.otm.indx = 1
            self.otm.gamma_bar = np.sqrt(11.0/12.0)*2*self.gam/(1+self.gam)
            self.otm.pf = np.pi*11.0**1.5*(1+2.0*self.gam**3)/(48*(1+self.gam)**3)
            self.otm.ratio_a = np.sqrt(12.0/11.0)*(1+self.gam)*0.5

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
