"""
:module: UnitTest_Pt3O4_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{Pt}_3\\mbox{O}_4` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.Pt3O4_lat as Pt
import unittest


class TestPt3O4(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.p_exp = 12
        self.mat_g = [0.12, 0.2, 0.3, 0.4, sqrt(2.0) - 1, 0.5, 0.6, 0.8, 0.9]
        self.a_nn = 3.0
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < sqrt(2.0) - 1:
                pf = np.pi * (3 + 4 * self.mat_g[i] ** 3) / 24.
            else:
                pf = np.pi * (3.0 + 4.0 * self.mat_g[i] ** 3) / (6. * sqrt(2.0) * (1 + self.mat_g[i]) ** 3)
            pt3o4 = Pt.LatPt3O4Base14(self.l_box, 1.0, self.mat_g[i])
            
            self.assertAlmostEqual(np.abs(pf - pt3o4.pf()), 0.0)
            
    def test_contacts(self):
        for i in range(len(self.mat_g)):
            pt3o4 = Pt.LatPt3O4Base14(self.l_box, self.a_nn, self.mat_g[i])
            
            if self.mat_g[i] < sqrt(2.0) - 1:
                vec1 = pt3o4.get_v([0, 0])
                vec2 = pt3o4.get_v([0, 1])
                vec = vec1 - vec2
                norm = self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)
                
                vec1 = pt3o4.get_v([0, 2])
                vec2 = pt3o4.get_v([0, 3])
                vec = vec1 - vec2
                norm = self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)
                
                vec1 = pt3o4.get_v([0, 4])
                vec2 = pt3o4.get_v([0, 5])
                vec = vec1 - vec2
                norm = self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)
                
            else:    
                vec1 = pt3o4.get_v([0, 0])
                vec2 = pt3o4.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)
                
                vec1 = pt3o4.get_v([0, 1])
                vec2 = pt3o4.get_v([1, 7])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)
                
                vec1 = pt3o4.get_v([0, 5])
                vec2 = pt3o4.get_v([1, 5])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                self.assertAlmostEqual(val, 1.0)


if __name__ == '__main__':
    unittest.main() 
