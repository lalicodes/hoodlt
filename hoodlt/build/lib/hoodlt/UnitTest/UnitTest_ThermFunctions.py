""":module:: UnitTest_Dmatrix_Pressure

   :platform: Unix, Windows
   :synopsis: Defines the unit test for the enthalpic pressure

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as Fcc
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.ThermFunctions as ThFs
import unittest


class TestThermFunctions(unittest.TestCase):

    def setUp(self):
        alat = 1.0
        cut_off_1 = 3.1
        lsize_1 = 6
        self.pval_1 = 6
        lat_1 = Fcc.LatFccBase4(lsize_1, alat)
        pot_1 = Ov.OverrpCutMixt(self.pval_1, np.array([1.0]), np.array([[1.0]]), cut_off_1)
        dmat_1 = ThFs.ThermoFunction(lat_1, pot_1)
        self.dm = dmat_1.d_mat.dmatrix_real()
        self.tm = dmat_1.d_mat.tmatrix_real()
        self.eps = 1e-8
        
        cut_off = 22.1
        lsize = 32
        self.pval = 14
        self.lat = Fcc.LatFccBase4(lsize, alat)
        self.pot = Ov.OverrpCutMixt(self.pval, np.array([1.0]), np.array([[1.0]]), cut_off)
        self.dmat = ThFs.ThermoFunction(self.lat, self.pot)
        
        # the value of this coefficient is from A. Travesset, J. Chem. Phys. 141, 164501 (2014)
        self.coeff_lit = 5.0346944
        
        self.mat_vib = self.dmat.dmat_vib(np.array([0, 0, 0]))

    def test_tmatrix(self):
        indx = np.where(np.abs(self.dm) > self.eps)
        indt = np.where(np.abs(self.tm[indx]/self.dm[indx] + (2+self.pval_1)) > self.eps)
        
        self.assertTrue(indt[0].size == 0)
    
    def test_entropy(self):
        coeff = 1 + np.log((2*np.pi)**(-1.5)*self.lat.g_l()**(-(1+0.5*self.pval))) - self.dmat.entropy()
        self.assertAlmostEqual(coeff, self.coeff_lit, 3)
        
    def test_vibrational_matrix(self):
        diagl = np.abs(np.diag(self.mat_vib) - self.mat_vib[0, 0])
        for val in diagl:
            self.assertTrue(val < 1e-14)
            
        
if __name__ == '__main__':
    unittest.main()
