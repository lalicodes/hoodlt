"""
:module: AlB2_OTM
:platform: Unix, Windows
:synopsis: Defines the class implementing the OTM :math:`\\mbox{AlB}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

from hoodlt.Lattices.Lat2.AlB2_lat import LatAlB2Base3
import hoodlt.OTM.OTMobs as ObT


class OTMLatAlB2Base3(LatAlB2Base3):
    """ Implementation of the :math:`\\mbox{AlB}_2` lattice
    three particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         """

        self.otm = ObT.OTMObs()
        super(OTMLatAlB2Base3, self).__init__(l_value, a_nn_e, gamma)

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
