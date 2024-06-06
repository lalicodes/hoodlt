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
import hoodlt.Lattices.Lat2.cub_AB13_lat as Cu
import unittest


class TestcubAB13(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 3
        self.a_nn = 1.0
        self.p_exp = 12

        self.g_crit_1 = (1+2.0*sqrt(2))/3.0-sqrt(3.0+4*sqrt(2))/3.0
        self.g_crit_2 = 1.0/(sqrt(5+2*sqrt(2))-1)

        self.mat_g = [0.12, 0.18, self.g_crit_1, 0.3, 0.4, 0.5, self.g_crit_2, 0.6, 0.8, 0.9]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < self.g_crit_1:
                pf = np.pi*(1 + 13 * self.mat_g[i]**3)/6.0
            elif (self.mat_g[i] >= self.g_crit_1) and (self.mat_g[i] < self.g_crit_2):
                a0 = sqrt(2)*(3*self.mat_g[i]**2-2*self.mat_g[i]-1)/(4*self.mat_g[i]
                                                                     - sqrt(-2*self.mat_g[i]**2+12*self.mat_g[i]+6))
                pf = np.pi*(1.0+13*self.mat_g[i]**3)/(6.0*a0**3)
            else:
                pf = np.pi*(13.0+1.0/self.mat_g[i]**3)/(6*(1+sqrt(2))**3)

            ab13 = Cu.LatCubAB13Base14(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertAlmostEqual(pf, ab13.pf(), 12)
            
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            ab13 = Cu.LatCubAB13Base14(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < self.g_crit_1:
                vec1 = ab13.get_v([0, 0])
                vec2 = ab13.get_v([0, 0]) + ab13.get_a(0)
                vec = vec1-vec2
                norm = a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0)
                self.assertTrue(np.abs(1 - val) < 1.e-12)
            
            elif self.mat_g[i] >= self.g_crit_1 and (self.mat_g[i] < self.g_crit_2):
                vec1 = ab13.get_v([0, 0])
                vec2 = ab13.get_v([1, 3])
                vec = vec1-vec2
                norm = 0.5*(1.0 + self.mat_g[i])*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0)
                self.assertTrue(np.abs(1 - val) < 1.e-12)
                
                vec1 = ab13.get_v([1, 1])
                vec2 = ab13.get_v([1, 0])
                vec = vec1-vec2
                norm = self.mat_g[i]*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                
                self.assertAlmostEqual(val, 1.0)
                self.assertTrue(np.abs(1 - val) < 1.e-12)
                
            else:
                vec1 = ab13.get_v([1, 0])
                vec2 = ab13.get_v([1, 1])
                vec = vec1-vec2
                norm = self.mat_g[i]*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                self.assertAlmostEqual(val, 1.0)
                self.assertTrue(np.abs(1 - val) < 1.e-12)
            
                vec1 = ab13.get_v([1, 2])+ab13.get_a(1)
                vec2 = ab13.get_v([1, 1])
                vec = vec1-vec2
                norm = self.mat_g[i]*a_nn
                val = sqrt(np.dot(vec, vec))/norm
                
                self.assertAlmostEqual(val, 1.0)
                self.assertTrue(np.abs(1 - val) < 1.e-12)

          
if __name__ == '__main__':
    unittest.main()
