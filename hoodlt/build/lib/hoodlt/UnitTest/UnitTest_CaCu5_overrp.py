"""
:module: UnitTest_CaCu5_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{CaCu}_5` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat2.CaCu5_lat as Ca
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCaCu5(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12
        self.mat_g = [0.12, 2.0 / np.sqrt(3.0) - 1, 0.2, 0.5, (1 + 2 * np.sqrt(19)) / 15.0, 0.7, 1 /
                      (4.0 / np.sqrt(3) - 1.0), 0.8, 0.9]

        self.mat_val = [13603.29576, 13603.29576, 9369.09208, 1625.35333, 902.08803,
                        401.21834, 171.39728, 98.18647, 23.89022]
        
        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            fac = 1 + 5 * self.mat_g[i] ** 3
            if self.mat_g[i] < 2.0 / np.sqrt(3.0) - 1:
                pf = np.pi * fac / (3.0 * np.sqrt(3.0))
            elif (self.mat_g[i] >= 2.0 / np.sqrt(3.0) - 1) and (self.mat_g[i] < (1 + 2 * np.sqrt(19)) / 15.0):
                pf = 4.0 * np.pi * fac / (9.0 * np.sqrt(3) * (1 + self.mat_g[i]) ** 2)
            elif (self.mat_g[i] >= (1 + 2 * np.sqrt(19)) / 15.0) and (self.mat_g[i] < 1 / (4.0 / np.sqrt(3) - 1.0)):
                pf = np.pi * 8.0 * fac / (
                 9 * np.sqrt(3) * ((1 + self.mat_g[i]) ** 2) * np.sqrt(15.0 * self.mat_g[i] ** 2
                                                                       - 2 * self.mat_g[i] - 1))
            else:
                pf = np.pi * fac / (24 * np.sqrt(2) * self.mat_g[i] ** 3)

            CaCu5 = Ca.LatCaCu5Base6(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertTrue(np.abs(CaCu5.pf() - pf) < 1.e-12)
            
    def test_density(self):
        for i in range(len(self.mat_g)):
            CaCu5 = Ca.LatCaCu5Base6(self.l_box, self.a_nn, self.mat_g[i])
            val = 2.0 * sum(CaCu5.num_pnts()) / (np.sqrt(3) * CaCu5.l_box()[0] * CaCu5.l_box()[1] * CaCu5.l_box()[2])
           
            self.assertTrue(np.abs(CaCu5.g_l() - val) < 1.e-12)

    def test_energy(self):
        for i in range(len(self.mat_g)):
            CaCu5 = Ca.LatCaCu5Base6(self.l_box, self.a_nn, self.mat_g[i])
            vp_CaCu5 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
            val = 2 * Dm.DMatrix(CaCu5, vp_CaCu5).energy()
            
            self.assertTrue(np.abs(val - self.mat_val[i]) < 1.e-5)
            
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            CaCu5 = Ca.LatCaCu5Base6(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < 2.0 / np.sqrt(3.0) - 1:
                vec1 = CaCu5.get_a(0)
                norm = a_nn
                val = np.sqrt(np.dot(vec1, vec1))/norm
            elif (self.mat_g[i] >= 2.0 / np.sqrt(3.0) - 1) and (self.mat_g[i] < (1 + 2 * np.sqrt(19)) / 15.0):
                vec1 = CaCu5.get_v([0, 0])
                vec2 = CaCu5.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * a_nn
                val = np.sqrt(np.dot(vec, vec))/norm
            elif (self.mat_g[i] >= (1 + 2 * np.sqrt(19)) / 15.0) and (self.mat_g[i] < 1 / (4.0 / np.sqrt(3) - 1.0)):
                vec1 = CaCu5.get_v([1, 0])
                vec2 = CaCu5.get_v([1, 2])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val = np.sqrt(np.dot(vec, vec))/norm
            else:
                vec1 = CaCu5.get_v([1, 2])
                vec2 = CaCu5.get_v([1, 3])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val = np.sqrt(np.dot(vec, vec))/norm
            
            self.assertAlmostEqual(val, 1.0)
            self.assertTrue(np.abs(val - 1.0) < 1.e-5)


if __name__ == '__main__':
    unittest.main()
