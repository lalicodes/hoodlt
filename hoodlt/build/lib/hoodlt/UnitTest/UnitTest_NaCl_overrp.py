"""
:module: UnitTest_NaCl_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{NaCl}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.NaCl_lat as Na
import hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt as Ova
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestNaCl(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 9
        self.a_nn = 1.0
        self.p_exp = 12
        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.mat_g = [0.2, 0.3, 0.4, sqrt(2) - 1, 0.5, 0.6, 1.0]
        self.mat_sc = 6.20209
        self.mat_fcc = [12.13181, 12.13181, 12.13181, 12.13181, 5.98423, 2.75843, 0.18953]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] < sqrt(2) - 1:
                pf = np.pi * (1 + self.mat_g[i] ** 3) / (3.0 * sqrt(2.0))
            else:
                pf = 2 * np.pi * (1 + self.mat_g[i] ** 3) / (3.0 * (self.mat_g[i] + 1) ** 3)

            nacl = Na.LatNaClBase8(self.l_box, self.a_nn, self.mat_g[i])
            
            self.assertAlmostEqual(nacl.pf(), pf)
            
    def test_sc_energy(self):
        nacl = Na.LatNaClBase8(self.l_box, self.a_nn, self.mat_g[-1])
        vp_nacl = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)
        val = 2 * Dm.DMatrix(nacl, vp_nacl).energy()
        
        self.assertAlmostEqual(val, self.mat_sc, 5)
        
    def test_fcc_energy(self):
        for i in range(len(self.mat_g)):
            nacl = Na.LatNaClBase8(self.l_box, self.a_nn, self.mat_g[i])
            vp_nacl = Ova.OverrpCutVanishMixt(self.p_exp, [0, 0], self.sigma_2, self.eps_2, self.cut_off)
            val = 2 * 2 * Dm.DMatrix(nacl, vp_nacl).energy()
            
            self.assertAlmostEqual(val, self.mat_fcc[i], 5)
            
         
if __name__ == '__main__':
    unittest.main()
