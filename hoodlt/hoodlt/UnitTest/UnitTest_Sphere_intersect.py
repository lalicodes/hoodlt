"""
:module: UnitTest_Sphere_Intersect.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the Sphere Intersect

.. moduleauthor:: Alex Travesset <trvsst@meslab.gov> , February 2020
"""

from __future__ import division

import unittest
import numpy as np

import hoodlt.GeneralPacking.sphere_intersect_cube as sc
import hoodlt.Groups.Rotation as rt


class TestLatticeNeighbors(unittest.TestCase):

    def test_intersect(self):

        # cube of edge 1 and sphere of diameter 1 with the same origin
        cube_c = np.array([0.0, 0.0, 0.0])
        cube_edge = 1.0
        cube_or = rt.Rotate(np.array([0.0, 0.0, 0.0]))

        sphere_c = np.array([0.0, 0.0, 0.0])
        sphere_d = 1.0

        val = sc.cube_sphere_intersect(cube_c, cube_edge, cube_or, sphere_c, sphere_d)
        # they intersect
        self.assertTrue(val)

        # cube of edge 1 and sphere of diameter 1-:math:`\\varepsilon` with both centers separated by 1
        cube_c = np.array([0.0, 0.0, 0.0])
        cube_edge = 1.0
        cube_or = rt.Rotate(np.array([0.0, 0.0, 0.0]))

        sphere_c = np.array([1.0, 0.0, 0.0])
        sphere_d = 1.0-1e-8

        val = sc.cube_sphere_intersect(cube_c, cube_edge, cube_or, sphere_c, sphere_d)
        # they do not intersect
        self.assertFalse(val)

        # cube of edge 1 and sphere of diameter 1 separated less than :math:`\\frac{\\sqrt{3}+1}{2}`
        dist = 0.5*(np.sqrt(3)+1)-1e-5
        cube_c = np.array([0.0, 0.0, 0.0])
        cube_edge = 1.0
        cube_or = rt.Rotate(np.array([0.0, 0.0, 0.0]))

        sphere_c = dist*np.array([1.0, 1.0, 1.0])/np.sqrt(3)
        sphere_d = 1.0

        val = sc.cube_sphere_intersect(cube_c, cube_edge, cube_or, sphere_c, sphere_d)
        # they do intersect
        self.assertTrue(val)

        # cube of edge 1 and sphere of diameter 1 separated more than :math:`\\frac{\\sqrt{3}+1}{2}`
        dist = 0.5 * (np.sqrt(3) + 1) + 1e-5
        cube_c = np.array([0.0, 0.0, 0.0])
        cube_edge = 1.0
        cube_or = rt.Rotate(np.array([0.0, 0.0, 0.0]))

        sphere_c = dist * np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        sphere_d = 1.0

        val = sc.cube_sphere_intersect(cube_c, cube_edge, cube_or, sphere_c, sphere_d)
        # they do not intersect
        self.assertFalse(val)

        # rotated cube with sphere along the x-axis
        cube_c = np.array([0.0, 0.0, 0.0])
        cube_edge = 1.0
        cube_or = rt.Rotate(np.array([np.pi/4, 0.0, 0.0]))

        sphere_c = np.array([1.0, 0.0, 0.0])
        sphere_d = 1.0
        val = sc.cube_sphere_intersect(cube_c, cube_edge, cube_or, sphere_c, sphere_d)
        # they do intersect
        self.assertTrue(val)


if __name__ == '__main__':
    unittest.main()
