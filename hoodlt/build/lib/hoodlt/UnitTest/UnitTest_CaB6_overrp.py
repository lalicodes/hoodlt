"""
:module: UnitTest_CaB6_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{CaB}_6` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:      
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat2.CaB6_lat as Ca
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCaB6(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_g = [0.2, 0.3, 0.4, 1.0 / (1.0 + np.sqrt(2)), 0.5, 0.6, 0.7]
        self.mat_val = [840323534.41442, 6477584.23973, 224615.52074, 168990.26318, 17656.85710, 1980.33392, 311.44179]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] > 1.0 / (1 + np.sqrt(2)):
                pf = np.pi * (1 + 6 * self.mat_g[i] ** 3) / (6.0 * (1 + np.sqrt(2)) ** 3 * self.mat_g[i] ** 3)
            else:
                pf = np.pi * (1 + 6 * self.mat_g[i] ** 3) / 6.0

            CaB6 = Ca.LatCaB6Base7(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertTrue(np.abs(CaB6.pf() - pf) < 1.e-12)

    def test_density(self):
        for i in range(len(self.mat_g)):
            CaB6 = Ca.LatCaB6Base7(self.l_box, self.a_nn, self.mat_g[i])
            val = ((7.0 * self.l_box ** 3) / (CaB6.l_box()[0] * CaB6.l_box()[1] * CaB6.l_box()[2]))
            
            self.assertEqual(CaB6.g_l() - val, 0.0)
            
    def test_energy(self):
        for i in range(len(self.mat_g)):
            CaB6 = Ca.LatCaB6Base7(self.l_box, self.a_nn, self.mat_g[i])
            vp_CaB6 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
            val = 2*Dm.DMatrix(CaB6, vp_CaB6).energy()
            
            self.assertTrue(np.abs(self.mat_val[i] - val) < 1.e-5)


if __name__ == '__main__':
    unittest.main()
