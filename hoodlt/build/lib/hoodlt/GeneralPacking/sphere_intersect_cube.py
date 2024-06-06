"""
:module: sphere_intersect_cube
:platform: Unix, Windows
:synopsis: Determines if a cube intersects with a sphere

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np


def cube_sphere_intersect(cube_c, cube_edge, cube_orientation, sphere_c, sphere_diam):
    """Determines if there is intersection between a cube and a sphere

    :param cube_c: center of the cube
    :param cube_edge: cube edge length
    :param cube_orientation: Rotate or RotateFromAxis object
    :param sphere_c: sphere center
    :param sphere_diam: sphere diameter
    """

    # cube edge vertices
    e1 = 0.5*cube_edge*np.ones(3)
    s1 = np.dot(cube_orientation.space_to_body, cube_c-sphere_c)

    return axis_align_cube_sphere_intersect(e1, s1, sphere_diam)


def axis_align_cube_sphere_intersect(e1, s1, sphere_diam):
    """Determines if a cube aligned with the 3D cartesian axis and its center at (0,0,0) intersects with a given sphere

    :param e1: corner of the cube in the positive quadrant
    :param s1: coordinates of sphere center
    :param sphere_diam: diameter of the sphere
    """

    e2 = -e1

    dist = (0.5*sphere_diam)**2

    for ind in range(3):
        if s1[ind] < e2[ind]:
            dist -= (s1[ind]-e2[ind])**2
        elif s1[ind] > e1[ind]:
            dist -= (s1[ind]-e1[ind])**2
        else:
            pass

    return dist > 0
