"""
:module: UnitTest_cub_fccAB13_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{cubfccAB}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""


import unittest
import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.cub_fccAB13_lat as Cf


class TestfccAB13(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 3
        self.a_nn = 1.0
        self.p_exp = 12

        self.mat_g = [0.12, 0.18, 1.0 - sqrt(2.0 / 3.0), 0.3, 0.4, (sqrt(10) + 1) / 9.0, 0.5, 0.6, 0.8, 0.9]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < 1-sqrt(2.0/3.0):
                pf = np.pi*(1+13*self.mat_g[i]**3)/(3.0*sqrt(2))
            elif (self.mat_g[i] >= 1-sqrt(2.0/3.0)) and (self.mat_g[i] < (sqrt(10)+1)/9.0):
                a0 = (4*self.mat_g[i]+2-6*(self.mat_g[i])**2)/(sqrt(2)*(-2*self.mat_g[i]+sqrt(-2*self.mat_g[i]**2+4*self.mat_g[i]+2)))
                pf = np.pi*(4+52*self.mat_g[i]**3)/(6.0*a0**3)
            else:
                pf = np.pi*(13.0+1.0/self.mat_g[i]**3)/(81.0*sqrt(2))

            fcc_ab13 = Cf.LatCubfccAB13Base56(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertAlmostEqual(fcc_ab13.pf(), pf, 12)
            
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            fcc_ab13 = Cf.LatCubfccAB13Base56(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < 1 - sqrt(2.0 / 3.0):
                vec1 = fcc_ab13.get_v([0, 0])
                vec2 = fcc_ab13.get_v([0, 1])
                vec = vec1 - vec2
                norm = a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                val2 = sqrt(np.dot(vec, vec)) / norm
            
            elif (self.mat_g[i] >= 1.0 - sqrt(2.0 / 3.0)) and (self.mat_g[i] < (sqrt(10) + 1) / 9.0):
                vec1 = fcc_ab13.get_v([0, 1])
                vec2 = fcc_ab13.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            
                vec1 = fcc_ab13.get_v([1, 26])
                vec2 = fcc_ab13.get_v([1, 0])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val2 = sqrt(np.dot(vec, vec)) / norm
            
            else:
                vec1 = fcc_ab13.get_v([1, 0])
                vec2 = fcc_ab13.get_v([1, 21])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            
                vec1 = fcc_ab13.get_v([1, 26])
                vec2 = fcc_ab13.get_v([1, 0])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val2 = sqrt(np.dot(vec, vec)) / norm
                
            self.assertAlmostEqual(np.abs(val - 1.0), 0.0, 12)
            self.assertAlmostEqual(np.abs(val2 - 1.0), 0.0, 12)
            
            
if __name__ == '__main__':
    unittest.main()
