"""
:module: ReciprocalVectors
:platform: Unix, Windows
:synopsis: calculates the reciprocal vectors of the unit cell and the maximum cut-off for minimal image convention

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2015
"""
from __future__ import division

import numpy as np
import numpy.linalg as la


def reciprocal_vectors(lat):
    """Returns the reciprocal vectors

    :param lat: DMixtLattice
    :return: numpy array with reciprocal vectors
    :rtype: numpy array
    
    Example
    
    .. code-block:: python
    
        print(reciprocal_vectors(lat))
        # returns a 3x3 numpy array containing the three reciprocal vectors
    """

    b = np.zeros((3, 3))

    # reciprocal space
    b[0] = np.cross(lat.get_a(1), lat.get_a(2)) / np.dot(lat.get_a(0), np.cross(lat.get_a(1), lat.get_a(2)))
    b[1] = np.cross(lat.get_a(2), lat.get_a(0)) / np.dot(lat.get_a(1), np.cross(lat.get_a(2), lat.get_a(0)))
    b[2] = np.cross(lat.get_a(0), lat.get_a(1)) / np.dot(lat.get_a(2), np.cross(lat.get_a(0), lat.get_a(1)))

    return b


def maximum_cutoff(lat):
    """Returns the maximum cut off for a cell that uses a restricted minimal distance convention

    :param lat: DMixtLattice
    :return: maximum cutoff
    :rtype: float
    
    Example
    
    .. code-block:: python
    
        print(maximum_cutoff(lat))
        # returns the maximum cutoff possible for the unit cell
    """

    b_vector = reciprocal_vectors(lat)

    values = lat.l/np.array([la.norm(b_vector[0]), la.norm(b_vector[1]), la.norm(b_vector[2])])

    return 0.5*np.amin(values)
