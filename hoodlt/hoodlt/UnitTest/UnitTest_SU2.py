"""
:module: UnitTest_SU2.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the SU2 class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> , May 2017
"""


import hoodlt.Groups.SU as SU
import numpy as np
import unittest


class TestSU2(unittest.TestCase):

    def test_quaternion_sum(self):
        """tests sums
        """

        ang1 = np.pi/3.0
        ang2 = np.pi/5.0
        vec = np.array([np.sin(ang2), np.cos(ang2), 0.0])
        quat1 = SU.SU2(ang1, vec).rot
        quat = quat1*quat1
        self.assertTrue(np.allclose(quat.norm(), 1.0))

if __name__ == '__main__':
    unittest.main()
