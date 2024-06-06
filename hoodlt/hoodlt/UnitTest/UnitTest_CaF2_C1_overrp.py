"""
:module: UnitTest_CaF2_C1_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{CaF}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  - rewrote test in proper unit test format
"""

import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.CaF2_lat as Ca
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt as Ova
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCaF2C1(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 9
        self.a_nn = 1.0
        self.p_exp = 12

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_g = [0.1, 0.2, 0.6, 1.0]
        self.mat_val_C1 = [2188.98061, 2188.98061, 88.58278, 6.08733]

        self.mat_val_fcc = [12.13181, 12.1318, 0.49092, 0.03372]
        self.mat_val_sc = [396.93740, 396.93740, 16.06307, 1.10382]
        
        self.vp_CaF2 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
        
        self.vp_CaF2_fcc = Ova.OverrpCutVanishMixt(self.p_exp, [0, 0], self.sigma_2, self.eps_2, self.cut_off)
        self.vp_CaF2_sc = Ova.OverrpCutVanishMixt(self.p_exp, [1, 1], self.sigma_2, self.eps_2, self.cut_off)
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < (sqrt(6.0) / 2 - 1):
                pf = np.pi * (1 + 2 * self.mat_g[i] ** 3) / (3 * sqrt(2))
            else:
                pf = np.pi * sqrt(3.0) * (1 + 2.0 * self.mat_g[i] ** 3) / (4.0 * (1 + self.mat_g[i]) ** 3)

            CaF2 = Ca.LatCaF2Base12(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertTrue(np.abs(pf - CaF2.pf()) < 1.e-12)
            
    def test_energy_C1(self):
        for i in range(len(self.mat_g)):
            CaF2 = Ca.LatCaF2Base12(self.l_box, self.a_nn, self.mat_g[i])
            val = 2 * Dm.DMatrix(CaF2, self.vp_CaF2).energy()
            
            self.assertTrue(np.abs(val - self.mat_val_C1[i]))

    def test_energy_fcc(self):
        for i in range(len(self.mat_g)):
            CaF2 = Ca.LatCaF2Base12(self.l_box, self.a_nn, self.mat_g[i])
            val = 3 * 2 * Dm.DMatrix(CaF2, self.vp_CaF2_fcc).energy()
            
            self.assertTrue(np.abs(val - self.mat_val_fcc[i]) < 1.e-5)
            
    def test_energy_sc(self):
        for i in range(len(self.mat_g)):
            CaF2 = Ca.LatCaF2Base12(self.l_box, self.a_nn, self.mat_g[i])
            val = 1.5 * 2 * Dm.DMatrix(CaF2, self.vp_CaF2_sc).energy()
            
            self.assertTrue(np.abs(val - self.mat_val_sc[i]) < 1.e-5)
            
            
if __name__ == '__main__':
    unittest.main()
