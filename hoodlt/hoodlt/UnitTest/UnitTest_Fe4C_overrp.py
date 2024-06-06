"""
:module: UnitTest_cub_AB13_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{cubAB}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.Fe4C_lat as Fe
import unittest


class TestFe4C(unittest.TestCase):
     
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12
        self.gamma_crit = (sqrt(3.0)-1)/(sqrt(3.0/2.0)+1)
        self.mat_g = [0.12, 0.2, self.gamma_crit, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < self.gamma_crit:
                pf = np.pi*(1+4*self.mat_g[i]**3)/6.0
            else:
                pf = 0.5*np.pi*sqrt(3)*(4.0+1.0/self.mat_g[i]**3)/(sqrt(3.0/2.0)+1.0+1.0/self.mat_g[i])**3

            fe4C = Fe.LatCFe4Base5(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertAlmostEqual(pf, fe4C.pf(), 12)
            
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            fe4C = Fe.LatCFe4Base5(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < self.gamma_crit:
                vec1 = fe4C.get_v([0, 0])
                vec2 = fe4C.get_a(0)
                vec = vec1 - vec2
                norm = a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0, 8)
            else:
                vec1 = fe4C.get_v([1, 0])
                vec2 = fe4C.get_v([1, 1])
                vec = vec1 - vec2
                norm = self.mat_g[i]*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0, 8)
                
                vec1 = fe4C.get_v([1, 0])
                vec2 = fe4C.get_v([0, 0])
                vec = vec1 - vec2
                norm = 0.5*(1 + self.mat_g[i])*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0, 8)
            
            
if __name__ == '__main__':
    unittest.main()
