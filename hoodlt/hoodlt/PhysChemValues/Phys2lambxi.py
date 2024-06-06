"""
:module: Phys2lambxi
:platform: Unix, Windows
:synopsis: Converts physical parameters into lambda and xi

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2016
"""

import numpy as np


def phys2lambxi(max_length, max_area, rad, ma):
    """
        Returns the :math:`\\lambda` and :math:`\\xi` parameters from :math:`L, R, \\sigma, A_0`

        :param max_length: max_length :math:`L`
        :param max_area: maximum area density :math:`A_0`
        :param rad: nanoparticle radius :math:`R`
        :param ma: grafting density (in chains/nm square) :math:`\\sigma`
        :return: parameters :math:`\\lambda`, :math:`\\xi`
        :rtype: numpy.ndarray
        """

    xi = ma / max_area

    lam = max_length / rad

    return np.array([lam, xi])
