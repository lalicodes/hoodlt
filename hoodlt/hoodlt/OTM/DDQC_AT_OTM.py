"""
:module: DDQC_AT_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`DDQC/OTM` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
"""

from hoodlt.Lattices.Lat2.DDQC_AT_lat import LatDDQCAT
import hoodlt.OTM.OTMobs as ObT


class OTMLatDDQC_AT_OTM(LatDDQCAT):
    """
        Implementation of the :math:`DDQC` quasi-lattice three particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatDDQC_AT_OTM, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        self.otm.indx = 0
        self.otm.gamma_bar = self.gam
        self.otm.pf = self.pf()

        return self.otm
