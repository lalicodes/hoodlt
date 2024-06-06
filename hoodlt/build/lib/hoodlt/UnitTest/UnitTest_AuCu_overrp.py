
"""
:module: UNitTest_AuCu_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the :math:`\\mbox{AuCu}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:       
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""

import numpy as np
import hoodlt.Lattices.Lat2.AuCu_lat as Au
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestAuCu(unittest.TestCase):

    def setUp(self):
        self.cut_off = 3.1
        self.l_box = 9
        self.a_nn = 1.0
        self.p_exp = 12
        self.l_lamb = 2.0
        
        self.sigma_2 = np.array([1.0, 1.0])
        self.eps_2 = np.array([[1.0, 1.0], [1.0, 1.0]])
        
        self.mat_g = [0.2, 0.4, 0.5, 0.7, np.sqrt(3)-1, 0.734, 0.74, 0.85, 1.0]
        self.mat_val = [51.20936, 51.20936, 51.20936, 51.20936, 51.20936, 50.52411, 48.49107, 24.77332, 12.13181]
        self.vp_AuCu = Ov.OverrpCutMixt(self.p_exp, self.sigma_2, self.eps_2, self.cut_off)

    def test_density(self):
        for i in range(len(self.mat_g)):
            AuCu = Au.LatAuCuBase4(self.l_box, self.a_nn, self.mat_g[i])
            ls = AuCu.l_box()
            rho = np.sum(AuCu.num_pnts())/(ls[0]*ls[1]*ls[2])
            
            self.assertTrue(np.abs(AuCu.g_l() - rho) < 1.e-10)
            
    def test_pf(self):
        for i in range(len(self.mat_g)):
            if self.mat_g[i] > np.sqrt(3)-1:
                pf = np.pi*(1 + self.mat_g[i]**3)/(3.0*2.0*np.sqrt(2.0)*np.sqrt(0.5*self.mat_g[i]**2+self.mat_g[i]-0.5))
            else:
                pf = np.pi*(1+self.mat_g[i]**3)/6.0
              
            AuCu = Au.LatAuCuBase4(self.l_box, self.a_nn, self.mat_g[i])
            self.assertTrue(np.abs(pf - AuCu.pf()) < 1.e-8)
            
    def test_energies(self):
        for i in range(len(self.mat_g)):
            AuCu = Au.LatAuCuBase4(self.l_box, self.a_nn, self.mat_g[i])
            val1 = self.mat_val[i]
            val2 = 2*Dm.DMatrix(AuCu, self.vp_AuCu).energy()
            
            self.assertTrue(np.abs(val1 - val2) < 1.e-5)
            
    def test_a_nn(self):
        for i in range(len(self.mat_g)):
            AuCu = Au.LatAuCuBase4(self.l_box, self.a_nn, self.mat_g[i])
            AuCu2 = Au.LatAuCuBase4(self.l_box, self.a_nn * self.l_lamb, self.mat_g[i])
            val2 = 2*Dm.DMatrix(AuCu, self.vp_AuCu).energy() 
            val3 = 2*Dm.DMatrix(AuCu2, self.vp_AuCu).energy()*self.l_lamb**self.p_exp
            
            self.assertTrue(np.abs(val2 - val3) < 0.2)
            
    def test_g_l(self):
        for i in range(len(self.mat_g)):
            AuCu = Au.LatAuCuBase4(self.l_box, self.a_nn, self.mat_g[i])
            AuCu2 = Au.LatAuCuBase4(self.l_box, self.a_nn * self.l_lamb, self.mat_g[i])
            
            self.assertTrue(np.abs(AuCu.pf() - AuCu2.pf()) < 1.e-9)
        
           
if __name__ == '__main__':
    unittest.main()
