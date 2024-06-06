""":module:: UnitTest_bccAB6_overrp.py
   :platform: Unix, Windows
   :synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{bccAB}_6` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  - updated test to proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat2.bccAB6_lat as Ab
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest

class TestbccAB6(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12
        self.mat_g = [0.12, 0.2, np.sqrt(5.0 / 3.0) - 1, 0.3, 0.4, 1 / (np.sqrt(10) - 1), 0.6, 0.8, 0.9]
        self.mat_val = [163043.20303, 163043.20303, 163043.20303, 149994.37354, 61639.45169, 36503.35062, 1605.38852, 50.85281, 12.37324]
        
        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < np.sqrt(5.0 / 3.0) - 1:
                pf = np.pi * np.sqrt(3.0) * (1 + 6 * self.mat_g[i] ** 3) / 8.0
            elif (self.mat_g[i] >= np.sqrt(5.0 / 3.0) - 1) and (self.mat_g[i] < 1 / (np.sqrt(10) - 1)):
                pf = 5.0 * np.pi * np.sqrt(5.0) * (1 + 6.0 * self.mat_g[i] ** 3) / (24.0 * (1 + self.mat_g[i]) ** 3)
            else:
                pf = np.pi * np.sqrt(2) * (6 + 1.0 / self.mat_g[i] ** 3) / 96.0

            bccb6 = Ab.LatBccAB6Base14(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertTrue(np.abs(bccb6.pf() - pf) < 1.e-12)
            
    def test_volume_fraction(self):
        for i in range(len(self.mat_g)):
            bccb6 = Ab.LatBccAB6Base14(self.l_box, self.a_nn, self.mat_g[i])
            val = np.sum(bccb6.num_pnts()) / (bccb6.l_box()[0] * bccb6.l_box()[1] * bccb6.l_box()[2])
            
            self.assertEqual(val, bccb6.g_l())
            
    def test_contact(self):
        for i in range(len(self.mat_g)):
            bccb6 = Ab.LatBccAB6Base14(self.l_box, self.a_nn, self.mat_g[i])
            a_nn = 3.0
        
            if self.mat_g[i] < np.sqrt(5.0 / 3.0) - 1:
                vec1 = bccb6.get_v([0, 0])
                vec2 = bccb6.get_v([0, 1])
                vec = vec1 - vec2
                norm = a_nn
                val = np.sqrt(np.dot(vec, vec) / norm)
            
            elif (self.mat_g[i] >= np.sqrt(5.0 / 3.0) - 1) and (self.mat_g[i] < 1 / (np.sqrt(10) - 1)):
                vec1 = bccb6.get_v([0, 0])
                vec2 = bccb6.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * a_nn
                val = np.sqrt(np.dot(vec, vec) / norm)
            
            else:
                vec1 = bccb6.get_v([1, 2])
                vec2 = bccb6.get_v([1, 1])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val = np.sqrt(np.dot(vec, vec) / norm)
            
            self.assertTrue(np.abs(val - 1.0) > 1.e-5)
            
    def test_energies(self):
        for i in range(len(self.mat_g)):
            bccb6 = Ab.LatBccAB6Base14(self.l_box, self.a_nn, self.mat_g[i])
            vp_bccB6 = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
            val = 2 * Dm.DMatrix(bccb6, vp_bccB6).energy()
            
            self.assertTrue(np.abs(val - self.mat_val[i]) < 1.e-5)
        

if __name__ == '__main__':
    unittest.main()
