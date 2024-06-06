"""
:module: CubAB13_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{cubAB}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np

from hoodlt.Lattices.Lat2.cub_AB13_lat import LatCubAB13Base14

import hoodlt.OTM.OTMobs as ObT


class OTMLatCubAB13Base14(LatCubAB13Base14):
    """ Implementation of the :math:`\\mbox{cubAB}_{13}` lattice
    fourteen particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatCubAB13Base14, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam < self.gamma_crit[2]:
            self.otm.indx = 0
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()
        else:
            self.otm.indx = 1
            self.otm.gamma_bar = self.gamma_crit[2]
            self.otm.pf = np.pi*(1+13*self.gam**3)/(6.0*self.gamma_crit[2]**3*(1+np.sqrt(2))**3)
            self.otm.ratio_b = self.otm.gamma_bar/self.gam
            
        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm
