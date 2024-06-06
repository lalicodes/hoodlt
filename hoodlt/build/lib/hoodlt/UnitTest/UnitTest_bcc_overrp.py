"""
:module: UnitTest_bcc_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the bcc lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
.. history:       
..                Alex Upah <alexupah@iastate.edu>  July 2022
..                  -rewrote test in proper unit test format
"""

import unittest
import numpy as np
import hoodlt.Lattices.Lat1.bcc_lat_Mixt as Bc
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov


class Testbcc(unittest.TestCase):

    def setUp(self):
        self.cut_off = 8.1
        self.num_pnts = 20
        self.mat_p = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        self.mat_bcc = [22.64, 14.76, 12.25, 11.05, 10.36, 9.89, 9.56, 9.31, 9.11, 8.95, 8.82, 8.70, 8.61]
        self.mat_err = [2.5, 0.2, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004, 0.004]
        
        self.sigma = np.array([1.0])
        self.eps_mat = np.array([[1.0]])
        
        self.bcc_p = Bc.LatBcc(self.num_pnts, 1.0)
        self.bcc_b = Bc.LatBccBase2(self.num_pnts, 1.0)
        
        self.g_p = self.bcc_p.g_l()
        self.g_b = self.bcc_b.g_l()

    def test_density(self):
        self.assertEqual(self.g_p, self.g_b)
       
    def test_energy(self):
        for i in range(len(self.mat_p)):
            vp = Ov.OverrpCutMixt(self.mat_p[i], self.sigma, self.eps_mat, self.cut_off)
            val_p = 2.0 * Dm.DMatrix(self.bcc_p, vp).energy()
            val_b = 2.0 * Dm.DMatrix(self.bcc_b, vp).energy()
            
            self.assertAlmostEqual(val_p, val_b)
            self.assertTrue(np.abs(val_p - self.mat_bcc[i]) < self.mat_err[i])
            
    def test_packing_fraction(self):
        
        self.assertEqual(self.bcc_p.pf(), self.bcc_b.pf())
        self.assertTrue(np.abs(self.bcc_b.pf()-np.sqrt(3.0)*np.pi/8) < 1.e-12)


if __name__ == '__main__':
    unittest.main()
