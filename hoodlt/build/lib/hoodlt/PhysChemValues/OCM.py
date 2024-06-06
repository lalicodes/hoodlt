"""
:module: OCM
:platform: Unix, Windows
:synopsis: Defines the function defining the OCM result

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np


def ocm_val(lam, xi):
    """
    Returns the ocm spacing between NCs

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :return: nanoparticle separation according to the ocm result (:math:`\\tau`)
    :rtype: double
    """
    
    return 0.5*(1+lam)*(-1+np.sqrt(1+8.0*(1+3*lam*xi)/(1+lam)**3))

def adj_ocm_val(lam):
    """
    Returns the adjusted ocm spacing based on fitting of simulation data

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :return: nanoparticle separation according to adjusted ocm fit (:math:`\\tau`)
    """

    return -0.296075932959 * np.exp(-2.79349329143 * lam) + 1.30497993089



def ocm_val_201(lam):
    """
    Returns the adjusted ocm spacing based on fitting of simulation data

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :return: nanoparticle separation according to adjusted ocm fit (:math:`\\tau`)
    """

    return 1 + 7/(10.15*2)


def ocm_val_140(lam):
    """
    Returns the adjusted ocm spacing based on fitting of simulation data

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :return: nanoparticle separation according to adjusted ocm fit (:math:`\\tau`)
    """

    return 1 + 5/(9.285*2)


def ocm_val_1289(lam):
    """
    Returns the adjusted ocm spacing based on fitting of simulation data

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :return: nanoparticle separation according to adjusted ocm fit (:math:`\\tau`)
    """

    return 1 + 9/(20.15*2)



def ocm_val_angle(lam, xi):
    """
    Returns the ocm spacing between NCs

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :return: nanoparticle cone angle according to the ocm result (:math:`\\tau`)
    :rtype: double
    """

    return np.rad2deg(np.arccos(0.5 * (-1 + np.sqrt(1 + 8.0 * (1 + 3 * lam * xi) / (1 + lam) ** 3))))


def ocm_mod_val(lam, rat, xi):
    """
    Returns the angle of the ocm cone

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param rat: maximum length of the cone contianing the maximum distance divided by R.
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :return: nanoparticle separation according to the ocm-variable length result (:math:`\\tau`)
    :rtype: double
    """

    return 0.5 * rat * (-1 + np.sqrt(1 + 8.0 * (1 + 3 * lam * xi) / rat ** 3))


def ocm_mod_val_angle(lam, rat, xi):
    """
    Returns the angle of the cone NCs

    :param lam: Defined as :math:`\\lambda=\\frac{L}{R}` the ratio of maximum hydrocarbon length to nanoparticle core radius
    :param rat: maximum length of the cone contianing the maximum distance divided by R.
    :param xi: Defined as :math:`\\xi=\\frac{A_0}{A}=\\sigma A_0` the ratio of smallest to the actual molecular area.
    :return: cone angle according to the ocm-variable length result (:math:`\\tau`)
    :rtype: double
    """

    return np.rad2deg(np.arccos(0.5 * (-1 + np.sqrt(1 + 8.0 * (1 + 3 * lam * xi) / rat ** 3))))
