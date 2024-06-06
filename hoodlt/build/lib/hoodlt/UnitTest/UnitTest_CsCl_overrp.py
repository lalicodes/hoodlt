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
import hoodlt.Lattices.Lat2.CsCl_lat as Cs
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestCsCl(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 9
        self.a_nn = 1.0
        self.p_exp = 12

        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_p = [4, 5, 6, 7, 8, 9]
        self.mat_g = [0.3, 0.4, 0.7, sqrt(3) - 1, 0.8, 1.0]
  
        self.mat_bcc_pow = [17.46696, 13.94097, 12.08171, 11.01364, 10.34499, 9.89192]

        self.mat_sc_4 = [369676.86344, 369676.86344, 369676.86344, 369676.86344, 232959.13110, 65794.53941]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            CsCl = Cs.LatCsClBase2(self.l_box, self.a_nn, self.mat_g[i])
            
            if self.mat_g[i] > sqrt(3) - 1:
                pf = 0.5 * np.pi * sqrt(3) * (1 + self.mat_g[i] ** 3) / (1 + self.mat_g[i]) ** 3
            else:
                pf = np.pi * (1 + self.mat_g[i] ** 3) / 6.0
                
            self.assertAlmostEqual(pf, CsCl.pf())
    
    def test_energy(self):
        CsCl = Cs.LatCsClBase2(self.l_box, self.a_nn, 1.0)
        for i in range(len(self.mat_p)):
            vp_CsCl = Ov.OverrpCutMixt(self.mat_p[i], self.sigma_2, self.eps_2, self.cut_off)
            val = 2 * Dm.DMatrix(CsCl, vp_CsCl).energy()
            
            self.assertTrue(np.abs(val - self.mat_bcc_pow[i]), 1.e-5)
            
            
if __name__ == '__main__':
    unittest.main()
