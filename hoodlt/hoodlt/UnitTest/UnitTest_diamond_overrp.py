"""
:module:: UnitTest_diamond_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the diamond lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.diamond_lat_Mixt as Dl
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestDiamond(unittest.TestCase):

    def setUp(self):
        self.mat_p = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        self.cut_off = [3.1, 6.1]
        self.cut_off2 = 6.1
        self.num_pnts = 9
        self.sigma = np.array([1.0])
        self.eps_mat = np.array([[1.0]])
        
        self.dia_p = [Dl.LatDiamondBase2(self.num_pnts, 1.0), Dl.LatDiamondBase2(14, 1.0)]
        self.g_p = self.dia_p[0].g_l()
        self.dia_b = [Dl.LatDiamondBase8(self.num_pnts, 1.0), Dl.LatDiamondBase8(14, 1.0)]
        self.g_b = self.dia_b[0].g_l()
        
    def test_density(self):
        
        self.assertEqual(self.g_p, self.g_b)
        
    def test_vectors(self):
        for i in range(len(self.mat_p)):
            for j in range(2):
                vp = Ov.OverrpCutMixt(self.mat_p[i], self.sigma, self.eps_mat, self.cut_off[j])
                val1 = 2*Dm.DMatrix(self.dia_p[j], vp).energy()
                val2 = 2*Dm.DMatrix(self.dia_b[j], vp).energy()
            
                self.assertAlmostEqual(val1, val2, 12)
                
    def test_pf(self):
        for i in range(2):
            self.assertEqual(self.dia_p[i].pf(), self.dia_b[i].pf())
            self.assertAlmostEqual(np.abs(self.dia_p[i].pf() - np.sqrt(3)*np.pi/16), 0.0, 12)
        
        
if __name__ == '__main__':
    unittest.main()
