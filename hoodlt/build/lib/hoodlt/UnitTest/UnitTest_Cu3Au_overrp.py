"""
:module: UnitTest_Cu3Au_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{CaF}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.Cu3Au_lat as Cu
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCu3Au(unittest.TestCase):

    def setUp(self):
        self.cut_off = 4.1
        self.l_box = 9
        self.a_nn = 1.0
        self.p_exp = 12.0

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_p = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        self.mat_fcc = [14.36231339, 13.34229219, 12.79853643, 12.49184244, 12.31109575, 12.20088779,
                        12.13187301, 12.08772471, 12.05899158, 12.04002397, 12.02735482]
        self.mat_g = [0.3, 0.4, sqrt(2) - 1, 0.5, 0.6]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            Cu3Au = Cu.LatCu3AuBase4(self.l_box, self.a_nn, self.mat_g[i])

            if self.mat_g[i] > sqrt(2) - 1:
                pf = np.pi * sqrt(2) * (1 / 3.0 + self.mat_g[i] ** 3) / (1 + self.mat_g[i]) ** 3
            else:
                pf = np.pi * (1 + 3.0 * self.mat_g[i] ** 3) / 6.0
                
            self.assertAlmostEqual(pf, Cu3Au.pf(), 12)
            
    def test_energy(self):
        Cu3Au = Cu.LatCu3AuBase4(self.l_box, self.a_nn, 1.0)
        for i in range(len(self.mat_p)):
            vp_Cu3Au = Ov.OverrpCutMixt(self.mat_p[i], self.sigma_2, self.eps_2, self.cut_off)
            val = 2 * Dm.DMatrix(Cu3Au, vp_Cu3Au).energy()     
            
            self.assertAlmostEqual(val, self.mat_fcc[i])
            
            
if __name__ == '__main__':
    unittest.main()
