"""
:module: UnitTest_Dmatrix_Mixture

:platform: Unix, Windows
:synopsis: Defines another unit test for the D-matrix

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test form
"""


import numpy as np
import hoodlt.Lattices.Lat1.sc_lat_Mixt as Sc
import hoodlt.Lattices.Lat1.bcc_lat_Mixt as Bcc
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt as OvVan
import hoodlt.Lattices.Lat2.CsCl_lat as CsCl
import unittest


class TestDmatrixMixture(unittest.TestCase):
    
    def setUp(self):
        cut_off = 2.71
        lsize = 5
        alat = 1.0
        sigma = np.array([1.0, 1.0])
        eps_mat = np.array([[1.0, 1.0], [1.0, 1.0]])

        pot_mixture = Ov.OverrpCutMixt(6, sigma, eps_mat, cut_off)
        pot_simple = Ov.OverrpCutMixt(6, np.array([sigma[0]]), np.array([[eps_mat[0, 0]]]), cut_off)
        bcc = Bcc.LatBccBase2(lsize, alat)
        cscl = CsCl.LatCsClBase2(lsize, alat, 1.0)

        dmat_bcc = Dm.DMatrix(bcc, pot_simple)
        dmat_cscl = Dm.DMatrix(cscl, pot_mixture)
        self.dmatrix_bcc = dmat_bcc.dmatrix_real()
        self.dmatrix_cscl = dmat_cscl.dmatrix_real()
        self.diff = np.abs(self.dmatrix_bcc - self.dmatrix_cscl) 
        self.diff_err = self.diff[np.where(self.diff > 1e-8)]
        
        self.dmatrix_bcc_k = dmat_bcc.dmatrix_k()
        self.dmatrix_cscl_k = dmat_cscl.dmatrix_k()
        self.diff_k = np.abs(self.dmatrix_bcc_k - self.dmatrix_cscl_k)
        self.diff_err_k = self.diff_k[np.where(self.diff_k > 1e-8)]
        
        lsize = 8
        cut_off = 3.56

        sc = Sc.LatSC(lsize, alat*2.0/np.sqrt(3))
        cscl = CsCl.LatCsClBase2(lsize, alat, 1.0)

        pot_simple = Ov.OverrpCutMixt(6, np.array([sigma[0]]), np.array([[eps_mat[0, 0]]]), cut_off)
        pot_mixture = OvVan.OverrpCutVanishMixt(6, [0, 0], sigma, eps_mat, cut_off)

        self.dmat_sc = Dm.DMatrix(sc, pot_simple)
        self.dmat_cscl = Dm.DMatrix(cscl, pot_mixture)

        self.egy_sc = self.dmat_sc.energy()
        self.egy_cscl = 2.0*self.dmat_cscl.energy()

    def test_shape(self):
        self.assertEqual(self.diff_err.shape[0], 0.0)
        
    def test_shape_k(self):
        self.assertEqual(self.diff_err_k.shape[0], 0.0)
        
    def test_energy(self):
        self.assertAlmostEqual(self.egy_sc, self.egy_cscl, 12)


if __name__ == '__main__':
    unittest.main()
