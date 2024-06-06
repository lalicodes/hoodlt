"""
:module: Polyhedra_Geometry
:platform: Unix, Windows
:synopsis: Utility to calculate different geometric properties of a polyhedra

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December2016
"""

import numpy as np
import numpy.linalg as la


def tetravol(a, b, c, d):
    """Computes the volume of a tetrahedra

    :param a: coordinate of first point
    :param b: coordinate of second point
    :param c: coordinate of third point
    :param d: coordinate of fourth point
    :return: volume of the tetrahedra
    :rtype: float
    """
    return np.abs(np.dot((a-d), np.cross((b-d), (c-d))))/6.0


def triarea(a, b, c):
    """Calculates the area of a triangle

    :param a: coordinate of first point
    :param b: coordinate of second point
    :param c: coordinate of third point
    :return: area of the triangle
    :rtype: float
    """
    return 0.5*la.norm((np.cross((a-b), (a-c))))


def polyarea(phd):
    """Calculates the area of a polygon whose vertices are defined by ver

    :param phd: :class:`hoodlt.Utils.GeomClasses.PolyHedra`
    :return: area
    :rtype: float
    """

    vert = phd.vertices
    fac = phd.faces

    area = 0.0
    bar = np.zeros([3])

    for ind_1, v in enumerate(fac):
        vt = vert[v]

        for ind_2 in range(3):
            bar = np.sum(vt, axis=0) / float(vt.shape[0])

        for j in range(vt.shape[0] - 1):
            area += triarea(vt[j], vt[j + 1], bar)
        area += triarea(vt[-1], vt[0], bar)

    return area


def polyvolume(phd):
    """Computes the volume of a polyhedra

    :param phd: :class:`hoodlt.Utils.GeomClasses.PolyHedra`
    :return: Volume
    :rtype: float
    """

    vert = phd.vertices
    fac = phd.faces

    # center of mass of the polyhedra
    pnto = np.sum(vert, axis=0)/float(vert.shape[0])

    vol = 0.0
    bar = np.zeros([3])

    for ind_1, v in enumerate(fac):
        vt = vert[v]

        for ind_2 in range(3):
            bar = np.sum(vt, axis=0) / float(vt.shape[0])

        for j in range(vt.shape[0] - 1):
            vol += tetravol(pnto, vt[j], vt[j + 1], bar)
        vol += tetravol(pnto, vt[-1], vt[0], bar)

    return vol


def isoperimetric(phd):
    """Computes the isoperimetric ratio of a polyhedra

    :param phd: :class:`hoodlt.Utils.GeomClasses.PolyHedra`
    :return: Volume
    :rtype: float
    """

    return 36.0*np.pi*polyvolume(phd)**2/polyarea(phd)**3
