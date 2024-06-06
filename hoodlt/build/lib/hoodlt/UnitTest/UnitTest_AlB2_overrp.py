"""
:module: UNitTest_AlB2_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{AlB}_2` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat2.AlB2_lat as AlB
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestAlB2(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.lin_box = 9
        self.a_nn = 1.0
        
        self.sigma = np.array([1.0, 1.0])
        self.eps_mat = np.array([[1.0, 1.0], [1.0, 1.0]])
        
        self.mat_g = [0.45, np.sqrt(7.0/3.0)-1, 0.5400, 0.57, 1/np.sqrt(3.0), 0.63, 0.68]
        self.mat_ex = [10.4612056968, 13.6889001011, 13.7987041804, 14.2090735560, 14.3443071999,
                       13.4897432855,  12.9103145438]
        
        self.mat_pf = [0.7147880997, 0.7821117807, 0.7802171790, 0.7788805001, 0.7792050851, 0.7366044326, 0.7099270686]

    def test_pf(self):
        for i in range(len(self.mat_g)):
            alb2 = AlB.LatAlB2Base3(self.lin_box, self.a_nn, self.mat_g[i])
            self.assertTrue(np.abs(alb2.pf() - self.mat_pf[i]) < 1.e-8)
            
    def test_density(self):
        for i in range(len(self.mat_g)):
            alb2 = AlB.LatAlB2Base3(self.lin_box, self.a_nn, self.mat_g[i])
            g_dens = np.sum(alb2.num_pnts())/(np.prod(alb2.l)*alb2.vol_unit_cell())
            
            self.assertTrue(np.abs(alb2.g_l() - g_dens) < 1.e-9)
            
    def test_energies(self):
        for i in range(len(self.mat_g)):
            alb2 = AlB.LatAlB2Base3(self.lin_box, self.a_nn, self.mat_g[i])
            
            self.sigma[1] = self.mat_g[i]
            vp = Ov.OverrpCutMixt(6, self.sigma, self.eps_mat, self.cut_off)
            val = 2.0 * Dm.DMatrix(alb2, vp).energy()
            
            self.assertTrue(np.abs(val - self.mat_ex[i]) < 1.e-10)
    

if __name__ == '__main__':
    unittest.main()
