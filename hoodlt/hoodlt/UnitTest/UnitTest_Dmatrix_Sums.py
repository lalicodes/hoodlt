"""
:module: UnitTest_Dmatrix_Sums

:platform: Unix, Windows
:synopsis: Defines the unit test for the sums of the D-matrix

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 202
..                  -rewrote test in proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat2.AlB2_lat as AlB2
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import numpy.linalg as la
import unittest


class TestDmatrixSums(unittest.TestCase):
    
    def setUp(self):
        self.num_pnts = 9
        sigma = np.array([1.0, 1.0])
        eps_mat = np.array([[1.0, 1.0], [1.0, 1.0]])
        cut_off = 2.71
        self.pot = Ov.OverrpCutMixt(6, sigma, eps_mat, cut_off)

        self.mat_g = [0.1, 0.4, 0.6, 0.8]
        
        d1_v1 = np.array([-693.1857082, -1852.0170879, -1852.0170879])
        d1_v2 = np.array([-693.1857082, -1852.0170879, -1852.0170879])
        d1_v3 = np.array([-479.66442147, -1344.63751021, -1344.63751021])
        d1_v4 = np.array([-133.27866212, -187.71984286, -187.71984286])
        self.mat_d1 = [d1_v1, d1_v2, d1_v3, d1_v4]
        
        d2_v1 = np.array([1632.56272126, 6771.30757768, 6771.30757768])
        d2_v2 = np.array([1632.56272126, 6771.30757768, 6771.30757768])
        d2_v3 = np.array([1118.94488878, 4934.84208571, 4934.84208571])
        d2_v4 = np.array([350.30562272, 621.85012281, 621.85012281])
        self.mat_d2 = [d2_v1, d2_v2, d2_v3, d2_v4]
        
    def test_linear_sums(self):
        for gam in self.mat_g:
            alb2 = AlB2.LatAlB2Base3(self.num_pnts, 1.0, gam)
            dm = Dm.DMatrix(alb2, self.pot)
            mat2 = dm.coeff_linear()
            
            self.assertTrue(np.amax(np.abs(mat2)) < 1.e-11)
            
    def test_d1(self):
        for ind, gam in enumerate(self.mat_g):
            alb2 = AlB2.LatAlB2Base3(self.num_pnts, 1.0, gam)
            dm = Dm.DMatrix(alb2, self.pot)
            dmat = dm.dmatrix_real()
            mat2 = dm.coeff_d1
            
            self.assertTrue(la.norm(self.mat_d1[ind] - mat2) < 1.e-8)
            
    def test_d2(self):
        for ind, gam in enumerate(self.mat_g):
            alb2 = AlB2.LatAlB2Base3(self.num_pnts, 1.0, gam)
            dm = Dm.DMatrix(alb2, self.pot)
            dmat = dm.dmatrix_real()
            mat2 = dm.coeff_d2
            
            self.assertTrue(la.norm(self.mat_d2[ind] - mat2[1, 1]) < 1.e-8)


if __name__ == '__main__':
    unittest.main()
