"""
:module: OXM
:platform: Unix, Windows
:synopsis: Defines the function defining the OCM result

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np


def oxm_val(lam, xi, theta):
    """
    Returns the NC separation for an arbitrary aperture angle

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :param theta: Overlap angle, so that it interpolates between :math:`\\theta = 0 \\mbox{ OPM and } \\theta = \\theta_{max} \\mbox{ OCM}`
    :return: nanoparticle separation according to the oxm result (:math:`\\tau`)
    :rtype: double
    """

    nc = np.cos(theta)

    return ((2*nc**2/(1+nc))*(1+3*xi*lam))**(1/3.0)
