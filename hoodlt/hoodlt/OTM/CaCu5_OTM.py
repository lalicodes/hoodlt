"""
:module: CaCu5_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{CaCu}_5` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.CaCu5_lat import LatCaCu5Base6

import hoodlt.OTM.OTMobs as ObT


class OTMLatCaCu5Base6(LatCaCu5Base6):
    """ Implementation of the :math:`\\mbox{CaCu}_5` lattice
    six particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        gamma_min_2 = (8+np.sqrt(19))/(10*np.sqrt(3))
        gamma_min_1 = (1+2*np.sqrt(19))/(np.sqrt(3)*(8+np.sqrt(19)))

        self.otm = ObT.OTMObs()
        self.otm.gamma_crit = np.array([0.0, gamma_min_1, gamma_min_2, 1.0])
        super(OTMLatCaCu5Base6, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam > self.gamma_crit[3]:
            self.otm.indx = 0
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()

        elif (self.gam > self.otm.gamma_crit[2]) and (self.gam <= self.gamma_crit[3]):
            self.otm.indx = 1
            g_b = 4.0*self.gam/np.sqrt(3)-1
            self.otm.gamma_bar = g_b
            self.otm.pf = np.pi*8.0*(1+5.0*self.gam**3)/(9*np.sqrt(3)*(1+g_b)**2*np.sqrt(15*g_b**2-2*g_b-1))
            self.otm.ratio_b = g_b/self.gam

        elif (self.gam > self.gamma_crit[2]) and (self.gam <= self.otm.gamma_crit[2]):
            g_c = self.gamma_crit[2]
            eta_crit = 8.0*np.pi*(1+5*g_c**3)/(9*np.sqrt(3)*(1+g_c)**2*np.sqrt(15*g_c**2-2*g_c-1))
            self.otm.indx = 1
            self.otm.gamma_bar = g_c
            self.otm.pf = eta_crit*(1+5*self.gam**3)/(1+5*g_c**3)
            self.otm.ratio_b = g_c/self.gam

        elif (self.gam > self.otm.gamma_crit[1]) and (self.gam <= self.gamma_crit[2]):
            g_c = self.gamma_crit[2]
            eta_crit = 4.0*np.pi*(1+5*g_c**3)/(9*np.sqrt(3)*(1+g_c)**2)
            self.otm.indx = 1
            self.otm.gamma_bar = g_c
            self.otm.pf = eta_crit*(5.0+1.0/self.gam**3)/(5.0+1.0/g_c**3)
            self.otm.ratio_b = self.gam/g_c

        elif (self.gam > self.gamma_crit[1]) and (self.gam <= self.otm.gamma_crit[1]):
            g_b = np.sqrt(3.0)*self.gam/(2-np.sqrt(3)*self.gam)
            self.otm.indx = 1
            self.otm.gamma_bar = g_b
            self.otm.pf = 4.0*np.pi*(5.0+1.0/self.gam**3)*g_b**3/(9.0*np.sqrt(3)*(1+g_b)**2)
            self.otm.ratio_a = self.gam/g_b

        else:
            self.otm.indx = -1
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
