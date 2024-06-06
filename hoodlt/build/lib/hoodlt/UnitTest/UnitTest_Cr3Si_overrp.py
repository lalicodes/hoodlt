"""
:module: UnitTest_Cr3Si_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{CaF}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""

import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.Cr3Si_lat as Cr
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCr3Si(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 8
        self.a_nn = 1.0
        self.p_exp = 12

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_g = [0.12, 0.2, sqrt(5.0 / 3.0) - 1, 0.3, 0.4, 1 / (sqrt(5) - 1), 0.9, 1.0]

    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < sqrt(5.0 / 3.0) - 1:
                pf = np.pi * sqrt(3.0) * (1 + 3 * self.mat_g[i] ** 3) / 8.0
            elif (self.mat_g[i] >= sqrt(5.0 / 3.0) - 1) and (self.mat_g[i] < 1 / (sqrt(5) - 1)):
                pf = 5.0 * np.pi * sqrt(5.0) * (1 + 3.0 * self.mat_g[i] ** 3) / (24.0 * (1 + self.mat_g[i]) ** 3)
            else:
                pf = np.pi * (3.0 + 1.0 / self.mat_g[i] ** 3) / 24.0

            Cr3Si = Cr.LatSiCr3Base8(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertTrue(np.abs(Cr3Si.pf() - pf) < 1.e-12)
            
    def test_energy(self):
        Cr3Si = Cr.LatSiCr3Base8(self.l_box, self.a_nn, self.mat_g[-1])
        vp_Cr3Si = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
        val = 2 * Dm.DMatrix(Cr3Si, vp_Cr3Si).energy()
        
        self.assertAlmostEqual(np.abs(val), 3.6168999)
         
    def test_contacts(self):
        a_nn = 3.0
        for i in range(len(self.mat_g)):
            Cr3Si = Cr.LatSiCr3Base8(self.l_box, a_nn, self.mat_g[i])

            if self.mat_g[i] < sqrt(5.0 / 3.0) - 1:
                vec1 = Cr3Si.get_v([0, 0])
                vec2 = Cr3Si.get_v([0, 1])
                vec = vec1 - vec2
                norm = a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            elif (self.mat_g[i] >= sqrt(5.0 / 3.0) - 1) and (self.mat_g[i] < 1 / (sqrt(5) - 1)):
                vec1 = Cr3Si.get_v([0, 0])
                vec2 = Cr3Si.get_v([1, 0])
                vec = vec1 - vec2
                norm = 0.5 * (1.0 + self.mat_g[i]) * a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            else:
                vec1 = Cr3Si.get_v([1, 0])
                vec2 = Cr3Si.get_v([1, 1])
                vec = vec1 - vec2
                norm = self.mat_g[i] * a_nn
                val = sqrt(np.dot(vec, vec)) / norm
            
            self.assertAlmostEqual(val, 1.0)
            self.assertTrue(np.abs(1 - val) < 1.e-12)


if __name__ == '__main__':
    unittest.main()
