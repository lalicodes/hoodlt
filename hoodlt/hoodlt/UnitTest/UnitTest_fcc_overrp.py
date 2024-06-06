"""
:module: UnitTest_fcc_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the bcc lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as Fc
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import unittest


class Testfcc(unittest.TestCase):
    
    def setUp(self):
        self.sigma = np.array([1.0])
        self.eps_mat = np.array([[1.0]])
        
        self.mat_p = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        self.mat_fcc = [25.34, 16.97, 14.45, 13.36, 12.80, 12.49, 12.31, 12.20, 12.13, 12.09, 12.06, 12.04, 12.03]
        self.mat_err = [2.4, 0.2, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004, 0.004]
        
        num_pnts = 20
        self.fcc_p = Fc.LatFcc(num_pnts, 1.0)
        self.fcc_b = Fc.LatFccBase4(num_pnts, 1.0)
        self.g_p = self.fcc_p.g_l()
        self.g_b = self.fcc_b.g_l()
        self.cut_off = 8.1
        
    def test_pf(self):
        self.assertAlmostEqual(self.fcc_b.pf(), np.sqrt(2)*np.pi/6.0, 12)
        self.assertAlmostEqual(self.fcc_p.pf(), np.sqrt(2)*np.pi/6.0, 12)
        
    def test_energy(self):
        for i in range(len(self.mat_p)):
            vp = Ov.OverrpCutMixt(self.mat_p[i], self.sigma, self.eps_mat, self.cut_off)
            val1 = 2*Dm.DMatrix(self.fcc_p, vp).energy()
            val2 = 2*Dm.DMatrix(self.fcc_b, vp).energy()
            
            self.assertAlmostEqual(val1, val2, 8)
            self.assertTrue(np.abs(val1 - self.mat_fcc[i]) < self.mat_err[i])
        
if __name__ == '__main__':
    unittest.main()
