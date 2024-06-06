"""
:module: UnitTest_MgCu2_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{MgCu}_{2}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  - rewrote test in proper unit test format
"""


import numpy as np
from numpy import sqrt
import hoodlt.Lattices.Lat2.MgCu2_lat as Mg
import hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt as Ova
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestMgCu2(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 1.2126
        self.l_box = 3
        self.p_exp = 6
        self.a_nn = 1.0
        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])
        
        self.mat_g = [0.4, 0.5, 0.6, 0.7, 0.8, sqrt(2.0 / 3.0), 0.85, 0.87, 0.90, 0.99]
        self.mat_diamond = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.142481, 2.733196, 2.230135, 1.258853]
        
    def test_pf(self):
        for i in range(len(self.mat_g)):
            mgcu2 = Mg.LatMgCu2Base24(self.l_box, self.a_nn, self.mat_g[i])
            if self.mat_g[i] > sqrt(2.0 / 3.0):
                pf = np.pi * sqrt(2.0) * (1 / self.mat_g[i] ** 3 + 2.0) / 24.0
            else:
                pf = sqrt(3.0) * np.pi * (1 + 2 * self.mat_g[i] ** 3) / 16.0
                
            self.assertAlmostEqual(pf, mgcu2.pf(), 12)
            
    def test_energy(self):
        for i in range(len(self.mat_g)):
            mgcu2 = Mg.LatMgCu2Base24(self.l_box, self.a_nn, self.mat_g[i])
            vp_mgcu2 = Ova.OverrpCutVanishMixt(self.p_exp, [0, 0], self.sigma_2, self.eps_2, self.cut_off)
            val = 3 * 2 * Dm.DMatrix(mgcu2, vp_mgcu2).energy()
            
            self.assertAlmostEqual(val, self.mat_diamond[i], 6)
            

if __name__ == '__main__':
    unittest.main()
