"""
:module: ZnO_lat
:platform: Unix, Windows
:synopsis: Defines the classes implementing the :math:`\\mbox{ZnO}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test form
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.ZnO_lat as Zn
import unittest


class TestZnO(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.p_exp = 12
        self.mat_g = [0.12, 0.2, sqrt(3.0 / 2.0) - 1, 0.3, 0.4, 0.6, 0.8, 0.9]
        
        self.a_nn = 3.0
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < sqrt(3.0 / 2.0) - 1:
                pf = np.pi * (1 + self.mat_g[i] ** 3) / (3.0 * sqrt(2.0))
            else:
                pf = np.pi * sqrt(3) * (1.0 + self.mat_g[i] ** 3) / (4.0 * (1.0 + self.mat_g[i]) ** 3)

            zno = Zn.LatZnOBase4(self.l_box, 1.0, self.mat_g[i])
        
            self.assertAlmostEqual(zno.pf(), pf)
        
    def test_contacts(self):
        for i in range(len(self.mat_g)):
            zno = Zn.LatZnOBase4(self.l_box, self.a_nn, self.mat_g[i])

            if self.mat_g[i] < sqrt(3.0 / 2.0) - 1:
                vec1 = zno.get_v([0, 0])
                vec2 = zno.get_v([0, 0]) + zno.get_a(0)
                vec = vec1 - vec2
                norm = self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            else:
                vec1 = zno.get_v([0, 0])
                vec2 = zno.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1 + self.mat_g[i]) * self.a_nn
                val = sqrt(np.dot(vec, vec)) / norm
                
            self.assertAlmostEqual(val, 1.0)
            
            
if __name__ == '__main__':
    unittest.main()
