"""
:module: Displacements
:platform: Unix, Windows
:synopsis: Utility to calculate displacement matrix

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2015
"""

import numpy as np
import numpy.linalg as la


def eps():
    """Defines a small number

    :return: returns a small number of :math:`10^{-8}`
    :rtype: float
    """
    return 1e-8


def vibrations_cell(d_mat, ncor, lat):
    """Given a dynamical matrix object and lattice object returns the quantity :math:`\\langle u_0 u_{\\vec n} \\rangle`

    :param d_mat: D-matrix object :class:`hoodlt.D_matrix.Dmatrix_Mixt_D.DMatrix`
    :param ncor: value of :math:`\\vec n` the correlator value at which to be computed.
    :param lat: lattice object :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
    :return: matrix vibration matrix
    :rtype: numpy
    """

    l_size = d_mat.Box_size()
    v_u = float(l_size[0] * l_size[1] * l_size[2])
    matd = d_mat.get(0, 0, 0)
    m = matd.shape[0]

    a_matrix = np.zeros((m, m), dtype=complex)

    b = np.zeros((3, 3))

    # reciprocal space
    b[0] = np.cross(lat.get_a(1), lat.get_a(2)) / np.dot(lat.get_a(0), np.cross(lat.get_a(1), lat.get_a(2)))
    b[1] = np.cross(lat.get_a(2), lat.get_a(0)) / np.dot(lat.get_a(1), np.cross(lat.get_a(2), lat.get_a(0)))
    b[2] = np.cross(lat.get_a(0), lat.get_a(1)) / np.dot(lat.get_a(2), np.cross(lat.get_a(0), lat.get_a(1)))

    vecn = ncor[0] * lat.get_a(0) + ncor[1] * lat.get_a(1) + ncor[2] * lat.get_a(2)

    for l1 in range(l_size[0]):
        for l2 in range(l_size[1]):
            for l3 in range(l_size[2]):
                veck = l1 * b[0] / float(l_size[0]) + l2 * b[1] / float(l_size[1]) + l3 * b[2] / float(l_size[2])
                phase = 2 * np.pi * np.dot(veck, vecn)
                if not ((l1 == 0) and (l2 == 0) and (l3 == 0)):
                    temp_matrix = np.exp(1j * phase) * la.inv(d_mat.get(l1, l2, l3))
                    a_matrix += temp_matrix / v_u
                else:
                    a_matrix += mode_zero(d_mat.get(0, 0, 0)) / v_u

    return a_matrix


def mode_zero(b):
    """Helper function to compute zero models

    :param b: square matrix
    :return: matrix with zero modes ommitted
    :rtype: numpy array
    """

    m3 = b.shape[0]

    a0 = np.zeros((m3, m3), dtype=complex)

    if m3 == 3:
        return a0

    else:
        w, v = la.eigh(b)
        temp = np.zeros((m3, m3), dtype=complex)
        u = np.conjugate(v)

        for il in range(m3):
            if np.abs(w[il]) > eps():
                for iter1 in range(m3):
                    for iter2 in range(m3):
                        temp[iter1, iter2] = v[iter1, il] * u[iter2, il]
                a0 += (1 / w[il]) * temp

        return a0
