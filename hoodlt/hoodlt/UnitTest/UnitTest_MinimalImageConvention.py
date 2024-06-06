"""
:module: UnitTest_MinimalImageConvention.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the minimal distance convention

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May 2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""


import numpy as np
import numpy.linalg as la
import hoodlt.D_matrix.Dmatrix_Mixt_distances as Dm
import hoodlt.Lattices.Lat1.hcp_lat_Mixt as Hc
import unittest


class TestMinimalImageConvention(unittest.TestCase):

    def setUp(self):
        num_l = 10
        self.hcp = Hc.LatHcpBase2(num_l, 1.0)
        self.d_d = Dm.MinConv(self.hcp, 1.0)
        self.val_list1 = [0.4, 0.45, 0.6, 0.6, 0.95, -0.6, -0.95]
         
    def test_min_dist1(self):
        vec1 = 0.4 * self.hcp.l[0] * self.hcp.get_a(0)
        val1 = la.norm(vec1)
        val2 = la.norm(self.d_d.min_dist(vec1))
        
        self.assertAlmostEqual(val1, val2, 5)
        
    def test_min_dist2(self):
        for i in range(2):
            vec1 = 0.6 * self.hcp.l[i] * self.hcp.get_a(i)
            val1 = la.norm(vec1 - self.hcp.l[i] * self.hcp.get_a(i))
            val2 = la.norm(self.d_d.min_dist(vec1))
            
            self.assertAlmostEqual(val1, val2, 5)
        
    def test_min_dist3(self):
        vec1 = 0.45 * self.hcp.l[1] * self.hcp.get_a(1)
        val1 = la.norm(vec1)
        val2 = la.norm(self.d_d.min_dist(vec1))
        
        self.assertAlmostEqual(val1, val2, 5)
        
    def test_min_dist4(self):
        vec1 = 0.95 * self.hcp.l[0] * self.hcp.get_a(0) + 0.95 * self.hcp.l[1] * self.hcp.get_a(1)
        val1 = la.norm(vec1 - self.hcp.l[0] * self.hcp.get_a(0) - self.hcp.l[1] * self.hcp.get_a(1))
        val2 = la.norm(self.d_d.min_dist(vec1))
        
        self.assertAlmostEqual(val1, val2, 5)
        
    def test_min_dist5(self):
        vec1 = -0.6 * self.hcp.l[0] * self.hcp.get_a(0)
        val1 = la.norm(vec1 + self.hcp.l[0] * self.hcp.get_a(0))
        val2 = la.norm(self.d_d.min_dist(vec1))
        
        self.assertAlmostEqual(val1, val2, 5)
        
    def test_min_dist6(self):
        vec1 = -0.95 * self.hcp.l[0] * self.hcp.get_a(0) - 0.95 * self.hcp.l[1] * self.hcp.get_a(1)
        val1 = la.norm(vec1 + self.hcp.l[0] * self.hcp.get_a(0) + self.hcp.l[1] * self.hcp.get_a(1))
        val2 = la.norm(self.d_d.min_dist(vec1))
        
        self.assertAlmostEqual(val1, val2, 5)
        
        
if __name__ == '__main__':
    unittest.main()
