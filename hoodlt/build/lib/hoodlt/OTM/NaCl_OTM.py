"""
:module: NaCl_OTM
:platform: Unix, Windows
:synopsis: Defines the class implementing the OTM :math:`\\mbox{NaCl}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

from hoodlt.Lattices.Lat2.NaCl_lat import LatNaClBase8

import hoodlt.OTM.OTMobs as ObT


class OTMLatNaClBase8(LatNaClBase8):
    """
    Implementation of the OTM :math:`\\mbox{NaCl}` lattice eight particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatNaClBase8, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        """

        self.otm.indx = 0
        self.otm.gamma_bar = self.gam
        self.otm.pf = self.pf()

        return self.otm
