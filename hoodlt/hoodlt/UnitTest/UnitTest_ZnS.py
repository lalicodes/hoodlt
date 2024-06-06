"""
:module: UnitTest_ZnS_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{ZnS}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.ZnS_lat as Zn
import unittest


class TestZnS(unittest.TestCase):

    def setUp(self):
        self.cut_off = 4.1
        self.l_box = 9
        self.a_nn = 1.0
        self.mat_g = [0.1, 0.2, sqrt(3 / 2.0) - 1, 0.5, 0.7, 1.0]

    def test_pf(self):
        for i in range(len(self.mat_g)):
            zns = Zn.LatZnSBase8(self.l_box, self.a_nn, self.mat_g[i])
            if self.mat_g[i] > sqrt(3 / 2.0) - 1:
                pf = sqrt(3.0) * np.pi * (1 + self.mat_g[i] ** 3) / (4.0 * (1 + self.mat_g[i]) ** 3)
            else:
                pf = np.pi / (3.0 * sqrt(2.0)) * (1 + self.mat_g[i] ** 3)
                
            self.assertAlmostEqual(zns.pf(), pf)
            
            
if __name__ == '__main__':
    unittest.main()
