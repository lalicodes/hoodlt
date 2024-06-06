"""
:module: UnitTest_A15_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the A15 lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
.. history:
..                Alex Upah <alexupah@iastate.edu> June 2022
..                  - Converted test to unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat1.A15_lat_Mixt as Af
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class Test_A15_Overrp(unittest.TestCase):

    def setUp(self):
        self.sigma = np.array([1.0])
        self.num_pnts = 9
        self.eps_mat = np.array([[1.0]])
        
        self.mat_p = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        
        self.a15_p = Af.a15_p = Af.LatA15Base8(11, 1.0)
        
    def test_a_nn_definition(self):
        cut_off = 2.1
        a_nn1 = 1.5

        vp1 = Ov.OverrpCutMixt(6, self.sigma, self.eps_mat, cut_off)
        
        a15_p1 = Af.LatA15Base8(11, 1)
        
        a15_p2 = Af.LatA15Base8(11, a_nn1)
        vp2 = Ov.OverrpCutMixt(6, self.sigma, self.eps_mat, a_nn1 * cut_off)
        
        a_nn2 = 2.0
        a15_p3 = Af.LatA15Base8(11, a_nn2)
        vp3 = Ov.OverrpCutMixt(6, self.sigma, self.eps_mat, a_nn2 * cut_off)
        
        en_p1 = Dm.DMatrix(a15_p1, vp1).energy()
        en_p2 = Dm.DMatrix(a15_p2, vp2).energy()
        en_p3 = Dm.DMatrix(a15_p3, vp3).energy()
        
        self.assertTrue(np.abs(en_p1-a_nn1**6*en_p2) < 1.e-3)
        self.assertEqual(np.abs(en_p1 - a_nn1 ** 6 * en_p2), 1.1546319456101628e-14)
        
        self.assertTrue(np.abs(en_p1-a_nn2**6*en_p3) < 1.e-3)
        self.assertEqual(np.abs(en_p1-a_nn2**6*en_p3), 0.0)
        
    def test_A15_vector_definition(self):
        cut_off = 1.21
        value_list = [1.776356839400250e-15, 2.664535259100375e-15, 2.664535259100375e-15, 2.664535259100375e-15,
                      3.330669073875469e-15, 3.330669073875469e-15, 3.774758283725532e-15, 3.774758283725532e-15,
                      3.774758283725532e-15, 3.774758283725532e-15, 3.996802888650563e-15, 4.440892098500626e-15,
                      4.662936703425657e-15]
       
        for i in range(len(self.mat_p)):
            vp = Ov.OverrpCutMixt(self.mat_p[i], self.sigma, self.eps_mat, cut_off)
            ex_uval = .75 + 3 * (2.0 / np.sqrt(5)) ** self.mat_p[i]
            ex_comp = Dm.DMatrix(self.a15_p, vp).energy()
           
            self.assertTrue(ex_uval - ex_comp < 1.e-5)
            self.assertAlmostEqual(ex_uval - ex_comp, value_list[i])
    
    def test_A15_pf(self): 
             
        self.assertTrue(self.a15_p.pf(), np.pi/6.0)


if __name__ == '__main__':
    unittest.main()
