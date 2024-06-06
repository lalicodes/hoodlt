"""
:module: UnitTest_Quaternion.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the Quaternion class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> , May 2017
.. history:
..                Alex Upah <alexupah@iastate.edu> June 2022
..                  - Updated function usage to fit with rotation
"""


import hoodlt.Groups.Quaternion as Quat
import numpy as np
import unittest


class TestQuaternion(unittest.TestCase):

    def test_quaternion_sum(self):
        """tests sums
        """
        quat1 = Quat.Quaternion(np.array([1.0, 1.0, 1.0, 1.0]))
        quat2 = Quat.Quaternion(np.array([2.0, 2.0, 2.0, 2.0]))
        quat3 = quat1 + quat2
        self.assertTrue(np.allclose(quat3.norm(), 6.0))

    def test_multiplication(self):
        """ Tests multiplication 
        """
        quat1 = Quat.Quaternion(np.array([1.0, 3.0, 2.0, 1.0]))
        quat2 = Quat.Quaternion(np.array([1.0, 2.0, 3.0, 4.0]))

        quat_prod = quat1*quat2
        quat3 = Quat.Quaternion(np.array([-15.0, 10.0, -5.0, 10.0]))

        self.assertTrue(quat3 == quat_prod)

    def test_conjugate(self):
        """tests conjugate
        """
        quat1 = Quat.Quaternion(np.array([-15.0, 10.0, -5.0, 10.0]))
        quat2 = Quat.Quaternion(np.array([-15.0, -10.0, 5.0, -10.0]))

        self.assertEqual(quat1.conjugate(), quat2)

    def test_inverse(self):
        """tests inverse
        """
        quat1 = Quat.Quaternion(np.array([-15.0, 10.0, -5.0, 10.0]))
        quat_inv = quat1.inverse()

        prod = quat1*quat_inv

        self.assertTrue(np.allclose(prod.quat[0], np.array([1.0, 0.0, 0.0, 0.0])))

    def test_narray(self):
        """Tests numpy array
        """
        np1 = np.zeros([2, 2, 4])
        np1[:, :, :] = [-15.0, 10.0, -5.0, 10.0]
        quat1 = Quat.Quaternion(np1)
        quat_inv = quat1.inverse()
        prod = quat1 * quat_inv

        for ind1 in range(2):
            for ind2 in range(2):
                self.assertTrue(np.allclose(prod.quat[ind1, ind2], np.array([1.0, 0.0, 0.0, 0.0])))

    def test_subtraction(self):
        """Tests subtraction
        """
        ang = np.pi/3.0
        quat1 = Quat.Quaternion(np.array([np.cos(ang), 0.0, np.sin(ang), 0.0]))

        quat2 = quat1 - quat1

        self.assertEqual(quat2, Quat.Quaternion(np.array(np.zeros(4))))

    def test_division(self):
        """Tests division
        """
        ang = np.pi/3.0
        quat1 = Quat.Quaternion(np.array([np.cos(ang), 0.0, np.sin(ang), 0.0]))

        quat2 = quat1/quat1

        self.assertEqual(quat2, Quat.Quaternion(np.array([1.0, 0.0, 0.0, 0.0])))

    def test_rotation_matrix(self):
        """Tests static method rotation_matrix()
        """
        quat1 = np.array([1, 0, 0, 0])
        rot_mat1 = np.eye(3)

        quat2 = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])
        rot_mat2 = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])

        quat3 = np.array([0, 1, 0, 0])
        rot_mat3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])

        self.assertTrue(np.allclose(rot_mat1, Quat.rotation_matrix(quat1)))
        self.assertTrue(np.allclose(rot_mat2, Quat.rotation_matrix(quat2)))
        self.assertTrue(np.allclose(rot_mat3, Quat.rotation_matrix(quat3)))


if __name__ == '__main__':
    unittest.main()
