"""
:module: UnitTest_Li3Bi_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{Li}_{3}\\mbox{Bi}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.MgZn2_lat as Mg
import hoodlt.Lattices.Lat2.MgCu2_lat as MgCu
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt as Ova
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestMgZn2(unittest.TestCase):

    def setUp(self):
        self.cut_off = 2.1
        self.l_box = 5
        self.p_exp = 6
        self.a_nn = 1.0

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_g = [0.4, 0.5, 0.6, 0.7, 0.8, sqrt(2.0 / 3.0), 0.85, 0.87, 0.90, 0.99]
        self.mat_diamond = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.142481, 2.733196, 2.230135, 1.258853]
        
    def test_energies(self):
        for i in range(len(self.mat_g)):
            self.sigma_2[1] = self.mat_g[i]
            mgzn2 = Mg.LatMgZn2Base12(self.l_box, self.a_nn, self.mat_g[i])
            vp_mgzn2 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
            
            mgcu2 = MgCu.LatMgCu2Base24(self.l_box, self.a_nn, self.mat_g[i])
            vp_mgcu2 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
            
            val1 = 2 * Dm.DMatrix(mgzn2, vp_mgzn2).energy()
            val2 = 2 * Dm.DMatrix(mgcu2, vp_mgcu2).energy()
            
            self.assertTrue(np.abs(val1 - val2) < .02)
    
    def test_pf(self):
        for i in range(len(self.mat_g)):
            self.sigma_2[1] = self.mat_g[i]
            mgzn2 = Mg.LatMgZn2Base12(self.l_box, self.a_nn, self.mat_g[i])
            
            if self.mat_g[i] > sqrt(2.0 / 3.0):
                pf = np.pi * sqrt(2.0) * (1 / self.mat_g[i] ** 3 + 2.0) / 24.0
            else:
                pf = np.sqrt(3.0) * np.pi * (1 + 2 * self.mat_g[i] ** 3) / 16.0
                
            self.assertAlmostEqual(pf, mgzn2.pf())
            
    def test_diamond_energy(self):
        cut_off = 1.2126
        for i in range(len(self.mat_g)):
            mgzn2 = Mg.LatMgZn2Base12(self.l_box, self.a_nn, self.mat_g[i])
            vp_mgzn2 = Ova.OverrpCutVanishMixt(self.p_exp, [0, 0], self.sigma_2, self.eps_2, cut_off)
            val = 3 * 2 * Dm.DMatrix(mgzn2, vp_mgzn2).energy()
            
            self.assertAlmostEqual(val, self.mat_diamond[i], 5)
            
            
if __name__ == '__main__':
    unittest.main()
