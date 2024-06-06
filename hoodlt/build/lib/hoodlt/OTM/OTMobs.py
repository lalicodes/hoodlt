"""
:module: OTMObs
:platform: Unix, Windows
:synopsis: Defines the class implementing the OTM objects

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np


class OTMObs(object):
    """
    Defines the OTM structure that is returned in every OTM lattice object
    
    :ivar gamma_crit: values of :math:`\\gamma_c` and :math:`\\gamma^c`
    :ivar indx: this quantity is -1 if the branch does not exist, 0 if returns OPM, 1 for OTM solution with charge
    :ivar pf: modified packing fraction
    :ivar gamma_bar: :math:`\\bar{\\gamma}`
    :ivar ratio_a: relative radius of A particles :math:`\\frac{\\bar{r}_A}{r_A}`
    :ivar ratio_b: relative radius of B particles :math:`\\frac{\\bar{r}_B}{r_B}`
    """

    def __init__(self):
        """
        The constructor
        """

        self.gamma_crit = np.array([0.0, 1.0])
        self.indx = -1
        # this quantity is -1 if the branch does not exist, 0 for the OPM result 
        # and 1 for an OTM solution with non-trivial charge
        self.pf = 0.0
        # packing fraction
        self.gamma_bar = 0.0
        # gamma bar
        self.ratio_a = 1.0
        self.ratio_b = 1.0
        # relative radius ratio
