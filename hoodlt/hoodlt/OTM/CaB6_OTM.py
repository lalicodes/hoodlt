"""
:module: CaB6_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{CaB}_6` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.CaB6_lat import LatCaB6Base7

import hoodlt.OTM.OTMobs as ObT


class OTMLatCaB6Base7(LatCaB6Base7):
    """ Implementation of the :math:`\\mbox{CaB}_6` lattice
    seven particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatCaB6Base7, self).__init__(l_value, a_nn_e, gamma)

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
        else:
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = np.pi*(1+6*self.gam**3)/6.0
            self.otm.ratio_b = self.otm.gamma_bar/self.gam
            
        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
