"""
:module: UnitTest_Li3Bi_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{Li}_{3}\\mbox{Bi}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  - rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.Li3Bi_lat as Li
import unittest


class TestLi3Bi(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12
        self.mat_g = [0.12, 0.2, 0.5*sqrt(6.0)-1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < 0.5*sqrt(6.0)-1:
                pf = np.pi*(1+3*self.mat_g[i]**3)/(3.0*sqrt(2.0))
            else:
                pf = sqrt(3.0)*np.pi*(1 + 3*self.mat_g[i]**3)/(4.0*(1 + self.mat_g[i])**3)

            li3bi = Li.LatLi3BiBase16(self.l_box, self.a_nn, self.mat_g[i])
            self.assertAlmostEqual(pf, li3bi.pf(), 12)
            
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            li3bi = Li.LatLi3BiBase16(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < 0.5*sqrt(6.0)-1:
                vec1 = li3bi.get_v([0, 0])
                vec2 = li3bi.get_v([0, 1])
                vec = vec1 - vec2
                norm = a_nn
                val = sqrt(np.dot(vec, vec))/norm
            
            else:
                vec1 = li3bi.get_v([0, 0])
                vec2 = li3bi.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5*(1.0 + self.mat_g[i])*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                
            self.assertAlmostEqual(val, 1.0, 12)
           
           
if __name__ == '__main__':
    unittest.main()
