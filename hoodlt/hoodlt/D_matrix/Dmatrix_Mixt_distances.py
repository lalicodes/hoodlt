"""
:module: D_matrix_Mixt_distances
:platform: Unix, Windows
:synopsis: Defines distances

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
"""

from __future__ import division
import numpy as np
import abc
import hoodlt.Lattices.reciprocal_vectors as rv


class MinConv(object):
    """Calculation of distances within the minimal image convention

    """

    def __init__(self, lat, cut_off):
        """The constructor

        :param lat: DMixtLattice object
        :param cut_off: cut-off distance
        """

        self.a_v = np.transpose(lat.a_vector)*lat.l
        self.b_v = np.transpose(rv.reciprocal_vectors(lat).T/lat.l)

        max_cut = rv.maximum_cutoff(lat)

        # make sure the cutoff does not exceed the maximum value
        if max_cut < cut_off:
            print(self.a_v)
            raise ValueError('Cutoff is too large, use a cutoff smaller than %1.5f' % max_cut)

    def min_dist(self, r):
        """computes the minimum distance convention

        :param r: numpy array of the coordinates
        :return: distance within the minimal image convention
        :rtype: numpy.ndarray
        """

        # convert to fractional coordinates
        s_f = np.tensordot(self.b_v, r, axes=[1, 0])
        # return minimal image
        s_f += -np.rint(np.around(s_f, 15))
        # back to x,y,z coordinates
        return np.tensordot(self.a_v, s_f, axes=[1, 0])


class MinConvAbs(object):
    """Calculation of distances within the minimal image convention

    """

    @abc.abstractmethod
    def __init__(self, l_box):
        """The constructor

        :param l_box: box size
        """

        self.l_value = l_box

    @abc.abstractmethod
    def min_dist(self, r):
        """computes the minimum distance convention

        :param r: numpy array of the coordinates
        :return: distance within the minimal image convention
        :rtype: numpy.ndarray
        """
        return r


class MinConvRect(MinConvAbs):
    """Implements the distance for a rectangular (cubic, tetragonal, orthorhombic) lattice
    """

    def __init__(self, l_box):
        """The constructor

        :param l_box: box size
        """

        super(MinConvRect, self).__init__(l_box)

    @classmethod
    def dist_name(cls, name):
        return name == 'Rect'

    def min_dist(self, r):
        """computes the minimum distance convention

        :param r: numpy array of the coordinates
        :return: distance within the minimal image convention
        :rtype: numpy.ndarray
        """
        return r[:] - np.rint(r[:] / self.l_value[:]) * self.l_value[:]


class MinConvHexa(MinConvAbs):
    """Implements the distance for a rectangular (cubic, tetragonal, orthorhombic)
    """

    def __init__(self, l_box):
        """The constructor

        :param l_box: box size
        """
        super(MinConvHexa, self).__init__(l_box)

    @classmethod
    def dist_name(cls, name):
        return name == 'Hexa'

    def min_dist(self, r):
        """computes the minimum distance convention

        :param r: numpy array of the coordinates
        :return: distance within the minimal image convention
        :rtype: numpy.ndarray
        """

        w = r[:]/self.l_value[:]

        z = self.l_value[1]/self.l_value[0]

        rx = w[0]-np.rint(w[0]+z*w[1]/np.sqrt(3.0))+0.5*z*np.rint(2.0*w[1]/np.sqrt(3.0))
        ry = w[1]-np.rint(2.0*w[1]/np.sqrt(3.0))*np.sqrt(3)/2.0
        rz = w[2]-np.rint(w[2])

        return np.array([rx, ry, rz])*self.l_value
