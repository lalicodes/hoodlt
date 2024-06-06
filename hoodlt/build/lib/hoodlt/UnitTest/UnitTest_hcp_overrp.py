"""
:module: UnitTest_hcp_overrp
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the hcp lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, September2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat1.hcp_lat_Mixt as Hc
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class Testhcp(unittest.TestCase):
    
    def setUp(self):
        self.cut_off = 3.1
        self.num_pnts = 9
        self.sigma = np.array([1.0])
        self.eps_mat = np.array([[1.0]])
        
        #these values are from Ashcroft and Mermin
        
        self.mat_p = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        self.mat_fcc = [13.36, 12.80, 12.49, 12.31, 12.20, 12.13, 12.09, 12.06, 12.04, 12.03]
        
    def test_energy(self):
        hcp = Hc.LatHcpBase2(self.num_pnts, 1.0)
        for i in range(len(self.mat_p)):
            vp = Ov.OverrpCutMixt(self.mat_p[i], self.sigma, self.eps_mat, self.cut_off)
            val = 2.0 * Dm.DMatrix(hcp, vp).energy()
            
            self.assertTrue(np.abs(val - self.mat_fcc[i]) < 0.1)
            
    def test_hcp_energies(self):
        #this value is calculated by Stillinger
        hcp_6 = 14.454897
        a_nn = [1.0, 1.25]
        cut_off = [22.05, 25.1]
        val_list = [7.e-4, 5.e-2]
        
        for i in range(2):
            hcp = Hc.LatHcpBase2(51, a_nn[i])
            vp = Ov.OverrpCutMixt(6, self.sigma, self.eps_mat, cut_off[i])
            val = 2.0*Dm.DMatrix(hcp, vp).energy()
            val = val*a_nn[i]**6 
            self.assertTrue(np.abs(hcp_6 - val) < val_list[i])
        
    def test_pf(self):
        hcp = Hc.LatHcpBase2(51, 1.25)
        self.assertEqual(hcp.pf(), np.sqrt(2)*np.pi/6.0)
        
        
if __name__ == '__main__':
    unittest.main()
