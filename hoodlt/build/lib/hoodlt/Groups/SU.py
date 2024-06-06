"""
:module: SU
:platform: Unix, Windows, OS
:synopsis: Defines element of SU(2) within a quaternion representation

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  May 2017
"""

from __future__ import division
import numpy as np
import numpy.linalg as la
import hoodlt.Groups.Quaternion as Quat


class SU2(object):
    """
    Defines an element of :math:`SU(2)`
    """

    def __init__(self, ang, vec_n):
        """The Constructor
        
        :param ang: rotation angle
        :param vec_n: normal vector (axis of rotation)
        """

        eps = 1e-10
        if np.abs(la.norm(vec_n)-1.0) > eps:
            vec_n = vec_n/la.norm(vec_n)

        vec = np.zeros(4)
        vec[0] = np.cos(ang)
        vec[1:] = np.sin(ang)*vec_n[:]

        self.rot = Quat.Quaternion(vec)
