"""
:module: UnitTest_Groups_Rotation.py
:platform: Unix, Windows
:synopsis: Defines the unit test for Rotation Groups

.. moduleauthor:: Alex Travesset <trvsst@meslab.gov> , February 2020
"""

from __future__ import division

import unittest
import numpy as np

import hoodlt.Groups.Rotation as rt


class TestLatticeRotations(unittest.TestCase):

    def test_rotation(self):

        # this unit test checks the formula
        # math:`{\\vec r}^{\prime}={\\vec n}({\\vec n}{\\vec r}+[{\\vec r}-{|\vec n}({\\vec n}{\\vec r}]\cos(\alpha)+
        #       ({\\vec r}\\times {\\vec n} \\sin(\\alpha)'

        def explicit_rot(rot):
            """ this is a rotation from space to body frame"""
            n = rot.axis_rotation
            a1 = np.outer(n, n)
            a2 = np.eye(3)
            a3 = np.array([[0, n[2], -n[1]], [-n[2], 0, n[0]], [n[1], -n[0], 0]])
            return a1 + (a2-a1)*np.cos(rot.angle_rotation) + a3*np.sin(rot.angle_rotation)

        eul = np.array([np.pi/10, np.pi/5, np.pi/12])
        rot1 = rt.Rotate(eul)
        mat_rot1 = rot1.space_to_body
        mat_rot1_exp = explicit_rot(rot1)
        # check that both matrices are the same
        self.assertTrue(np.allclose(np.dot(mat_rot1, rot1.body_to_space), np.eye(3)))
        self.assertTrue(np.allclose(mat_rot1, mat_rot1_exp))

        eul = np.array([np.pi / 10, -np.pi / 5, np.pi / 12])
        rot1 = rt.Rotate(eul)
        mat_rot1 = rot1.space_to_body
        mat_rot1_exp = explicit_rot(rot1)
        # check that both matrices are the same
        self.assertTrue(np.allclose(np.dot(mat_rot1, rot1.body_to_space), np.eye(3)))
        self.assertTrue(np.allclose(mat_rot1, mat_rot1_exp))

        eul = np.array([np.pi / 10, -np.pi / 5, -np.pi / 12])
        rot1 = rt.Rotate(eul)
        mat_rot1 = rot1.space_to_body
        mat_rot1_exp = explicit_rot(rot1)
        # check that both matrices are the same
        self.assertTrue(np.allclose(np.dot(mat_rot1, rot1.body_to_space), np.eye(3)))
        self.assertTrue(np.allclose(mat_rot1, mat_rot1_exp))


    def test_rotation_from_axis(self):
        """
        Test the rotation from axis
        """

        naxis = np.array([1.0, 0.0, 0.0])
        ang = np.pi/12

        # define rotation
        rot1 = rt.RotateFromAxis(naxis, ang)

        # this rotation is the same as Euler angle \theta = angle
        eul = np.array([0.0, ang, 0.0])
        rot2 = rt.Rotate(eul)

        # test
        self.assertTrue(np.allclose(rot1.space_to_body, rot2.space_to_body))

        naxis = np.array([0.0, 0.0, 1.0])
        ang = np.pi / 4
        rot2 = rt.RotateFromAxis(naxis, ang)
        # this is a rotation of 45 degrees along the z-axis. the body coordinates of the vector
        vec_body = np.array([1, 0, 0])
        # should transform into the space coordinates
        vec_space = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0.0])

        self.assertTrue(np.allclose(np.dot(rot2.body_to_space, vec_body), vec_space))


if __name__ == '__main__':
    unittest.main()
