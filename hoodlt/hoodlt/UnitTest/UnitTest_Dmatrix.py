"""
:module: UnitTest_Dmatrix

:platform: Unix, Windows
:synopsis: Defines the unit test for the D-matrix

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test form
"""


import numpy as np
import hoodlt.Lattices.Lat2.AlB2_lat as AlB2
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import hoodlt.D_matrix.Dmatrix_Mixt_distances as Ds
import unittest


class TestDmatrix(unittest.TestCase):

    def setUp(self):
        self.num_pnts = 8
        self.cut_off = 2.71
        self.gam = 0.4
        self.p_exp = 6

        self.sigma = np.array([1.0, self.gam])
        self.eps_mat = np.array([[1.0, 1.0], [1.0, 1.0]])

        self.pot = Ov.OverrpCutMixt(self.p_exp, self.sigma, self.eps_mat, self.cut_off)
        self.lat = AlB2.LatAlB2Base3(self.num_pnts, 1.0, self.gam)

        self.dmat = Dm.DMatrix(self.lat, self.pot)
        
        self.rval = self.dmat.r_vec
        self.d_matrix = self.dmat.dmatrix_real()
        self.d_energy = self.dmat.energy()
        self.dist = Ds.MinConv(self.lat, self.pot.cut)
        
        self.delta_1 = np.array([0.0, 0.000001, 0.00001, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0008, 0.001])
        self.delta_2 = -np.array([0.0, 0.000001, 0.00001, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0008, 0.001])
        
    def test_sumrule(self):
        self.assertAlmostEqual(np.abs(np.amax(np.sum(self.dmat.dmatrix_real(), axis=(1, 2, 3, 4)))), 0.0, 12)
        
    def test_energies(self):
        ind1 = 0
        ind2 = 2
        typ1 = 0
        typ2 = 0

        if ind1 > self.lat.typ[0]-1:
            typ1 = 1

        if ind2 > self.lat.typ[0]-1:
            typ2 = 1
            
        for self.delt1, self.delt2 in zip(self.delta_1, self.delta_2):
            e0_exact = 0.0
            s_1 = np.copy(self.rval[:, ind1, :, :, :])
            r_1 = np.copy(s_1)
            r_1[0, :, :, :, :] += self.delt1
            r_1[0, ind1, 0, 0, 0] -= self.delt1
            r_1[0, ind2, 0, 0, 0] -= self.delt2
            s_2 = np.copy(self.rval[:, ind2, :, :, :])
            r_2 = np.copy(s_2)
            r_2[0, :, :, :, :] += self.delt2
            r_2[0, ind2, 0, 0, 0] -= self.delt2
            r_2[0, ind1, 0, 0, 0] -= self.delt1
            
            l1 = typ1
            
            for l2 in range(2):
                c_ind = np.arange(self.dmat.ncoord[l2], self.dmat.lcoord[l2])
                b_mat = r_1[:, c_ind, :, :, :]
                e0_exact += np.sum(self.pot.e1(self.dist.min_dist(b_mat), [l1, l2]))
                b_mat = s_1[:, c_ind, :, :, :]
                e0_exact += -np.sum(self.pot.e1(self.dist.min_dist(b_mat), [l1, l2]))

            l2 = typ2
            for l1 in range(2):
                c_ind = np.arange(self.dmat.ncoord[l1], self.dmat.lcoord[l1])
                b_mat = r_2[:, c_ind, :, :, :]
                e0_exact += np.sum(self.pot.e1(self.dist.min_dist(b_mat), [l1, l2]))
                b_mat = s_2[:, c_ind, :, :, :]
                e0_exact += -np.sum(self.pot.e1(self.dist.min_dist(b_mat), [l1, l2]))

            l1 = typ1
            l2 = typ2
            b_mat = r_1[:, ind2, 0, 0, 0]
            e0_exact -= self.pot.e1(self.dist.min_dist(b_mat), [l1, l2])
            b_mat = s_1[:, ind2, 0, 0, 0]
            e0_exact += self.pot.e1(self.dist.min_dist(b_mat), [l1, l2])

            e0_quad = self.d_matrix[3*ind1, 3*ind1, 0, 0, 0]*self.delt1**2
            e0_quad += self.d_matrix[3*ind2, 3*ind2, 0, 0, 0]*self.delt2**2 + 2 * self.delt2 * self.delt1 * \
                self.d_matrix[3*ind1, 3*ind2, 0, 0, 0]
            e0_quad *= 0.5
            
            self.assertAlmostEqual(e0_quad, e0_exact, 8)
            
            
if __name__ == '__main__':
    unittest.main()
