"""
:module: CFe4_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{Fe}_4\\mbox{C}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.Fe4C_lat import LatCFe4Base5
import hoodlt.OTM.OTMobs as ObT


class OTMLatCFe4Base5(LatCFe4Base5):
    """ Implementation of the :math:`\\mbox{Fe}_4\\mbox{C}` lattice
    five particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatCFe4Base5, self).__init__(l_value, a_nn_e, gamma)

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
        
        else:
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[1]
            self.otm.pf = np.pi*(1+4*self.gam**3)/6.0
            self.otm.ratio_b = self.otm.gamma_bar/self.gam
            
        if self.otm.pf > 1:
            self.otm.indx = -1
            
        return self.otm
