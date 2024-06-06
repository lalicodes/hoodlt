"""
:module: Lattice Coordinates
:platform: Unix, Windows
:synopsis: Utility to calculate conversion from lattice coordinates in cartesian to reduced and back

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December2016
"""

import numpy as np

from hoodlt.Lattices import reciprocal_vectors as rv


class LatCords(object):
    """
    Class to convert from cartesian to fractional coordinates
    """

    def __init__(self, lat, eps=1e-12):
        """The constructor

        :param lat: :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
        :param eps: precision value
        """

        self.prim = np.zeros([3, 3])
        self.inv = rv.reciprocal_vectors(lat)

        for ind in range(3):
            self.prim[ind] = lat.get_a(ind)

        self.eps = eps

    def fractional(self, v):
        """converts cartesian coordinates within eps into fractional coordinates

        :param v: vector in cartesian coordinates
        :return: vector in fractional coordinates
        :rtype: numpy array
        """

        val = np.dot(self.inv, v)

        return np.where(np.abs(val-np.rint(val)) > self.eps, val, np.rint(val))

    def cartesian(self, v):
        """converts fractional candidates into cartesian

        :param v: vector in fractional coordinates
        :return: vector in cartesian coordinates
        :rtype: numpy array
        """

        return np.dot(v, self.prim)

    def fractional_unitcell(self, v):
        """converts cartesian coordinates within eps into fractional coordinates within the unit cell

        :param v: vector in cartesian coordinates
        :return: vector in fractional coordinates
        :rtype: numpy array
        """

        # make sure the vector is within the unit cell
        vec = self.fractional(v) - np.trunc(self.fractional(v))

        # fix those values that are negative
        return np.where(vec < 0, 1+vec, vec)
