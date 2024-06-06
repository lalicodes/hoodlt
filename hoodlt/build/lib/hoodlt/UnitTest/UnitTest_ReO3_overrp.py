"""
:module: UnitTest_ReO3_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{ReO}_3` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.ReO3_lat as Re
import unittest


class TestReO3(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.p_exp = 12
        self.mat_g = [0.12, 0.2, sqrt(2.0)-1, 0.5, 0.6, 0.8, 0.9]
        self.a_nn = 3.0
        self.val_list = [0, 1, 2]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            pf = np.pi*(1 + 3.0*self.mat_g[i]**3) / (6.0*(1 + self.mat_g[i])**3)
            reo3 = Re.LatReO3Base4(self.l_box, 1.0, self.mat_g[i])
            
            self.assertAlmostEqual(reo3.pf(), pf)
            
    def test_contact(self):
        for i in range(len(self.mat_g)):
            reo3 = Re.LatReO3Base4(self.l_box, self.a_nn, self.mat_g[i])
            for j in self.val_list:
                vec1 = reo3.get_v([1, j])
                vec2 = reo3.get_v([0, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                
                self.assertAlmostEqual(val, 1.0)
                
                
if __name__ == '__main__':
    unittest.main() 
