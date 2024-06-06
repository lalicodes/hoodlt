"""
:module: OPM
:platform: Unix, Windows
:synopsis: Defines the function defining the OPM result

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2016
"""

import numpy as np


def opm_val(lam, xi):
    """
    Returns the opm spacing between NCs

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :return: nanoparticle separation according to the opm result (:math:`\\tau`)
    :rtype: double
        """

    return np.power(1+3*lam*xi, 1.0/3.0)
