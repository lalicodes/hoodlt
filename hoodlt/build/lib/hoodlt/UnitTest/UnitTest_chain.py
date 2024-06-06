"""
:module: UnitTest_Chain.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the chain_lat_Mixt class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> November 2018
"""

import unittest
import hoodlt.Lattices.Lat_dim1_1.chain_lat_Mixt as Chain
import numpy as np


class TestLatChain(unittest.TestCase):

    def setUp(self):
        self.l_val = 6
        self.a_nn = 1
        self.config = Chain.LatChain(self.l_val, self.a_nn)

    def test_l(self):
        self.assertEqual(self.config.l[0], self.l_val)
        self.assertEqual(self.config.l[1], 1)
        self.assertEqual(self.config.l[2], 1)
        self.assertEqual(((self.config.l_box()-np.array([self.l_val, 1, 1])*self.a_nn) == 0).all(), True)

    def test_a_vector(self):
        self.assertEqual(((self.config.a_vector[0]-np.array([1.0, 0.0, 0.0])*self.a_nn) == 0).all(), True)
        self.assertEqual(((self.config.a_vector[1]-np.array([0.0, 1.0, 0.0])*self.a_nn) == 0).all(), True)
        self.assertEqual(((self.config.a_vector[2]-np.array([0.0, 0.0, 1.0])*self.a_nn) == 0).all(), True)

    def test_type(self):
        self.assertEqual(self.config.typ, [1])

    def test_v_vector(self):
        self.assertEqual(((self.config.v_vector[0]-np.array([0.0, 0.0, 0.0])) == 0).all(), True)

    def test_num_pnts(self):
        self.assertEqual(self.config.num_pnts(), int(np.prod(self.l_val)))

    def test_gl(self):
        self.assertEqual(self.config.g_l(), 1)


if __name__ == '__main__':
    unittest.main()
