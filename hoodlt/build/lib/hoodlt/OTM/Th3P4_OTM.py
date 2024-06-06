"""
:module: Th3P4_OTM
:platform: Unix, Windows
:synopsis: Defines the class implementing the OTM :math:`\\mbox{Th}_3\\mbox{P}_4` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December2019
"""

from hoodlt.Lattices.Lat2.Th3P4_lat import LatTh3P4Base28
import hoodlt.OTM.OTMobs as ObT


class OTMLatTh3P4Base28(LatTh3P4Base28):
    """ Implementation of the :math:`\\mbox{Th}_3\\mbox{P}_4` lattice
    three particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatTh3P4Base28, self).__init__(l_value, a_nn_e, gamma)

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
