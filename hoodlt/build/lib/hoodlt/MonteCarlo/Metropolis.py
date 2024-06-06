"""
:module: Metropolis
:platform: Unix, Windows
:synopsis: Implements Monte Carlo Metropolis test

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
"""

from __future__ import division
import numpy as np


def metropolis_test(deltae, rnd):
    """Implements Metropolis test for acceptance/rejection
    
    :param deltae: difference in energy between two configurations (in units of :math:`k_B T`)
    :param rnd: random number uniformly distributed between [0, 1]
    :return: True if configurations swap, False otherwise
    :rtype: bool
    """

    swap = False

    if rnd < np.exp(-deltae):
        swap = True

    return swap
