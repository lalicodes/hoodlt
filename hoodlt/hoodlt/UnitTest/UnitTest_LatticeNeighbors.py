"""
:module: UnitTest_LatticeNeighbors.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the NanoInitGSD function

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> , April 2017
"""

from __future__ import division

import unittest
import numpy as np
import numpy.linalg as la

import hoodlt.Utils.LatticeNeighbors as ln
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as fcc
import hoodlt.Lattices.Lat2.MgZn2_lat as mgzn


class TestLatticeNeighbors(unittest.TestCase):

    def test_fcc(self):

        # the first distances of nearest neighbors are
        fcc_d = np.array([1.0, np.sqrt(2), np.sqrt(3), 2.0])
        # with multiplicity
        fcc_mult = np.array([12, 6, 24, 12])

        l_lat = 3
        # list of lattice constants
        a_lat = np.array([1/np.sqrt(2), 1.0, 3.0])
        # list of points
        lst = [0, 1, -1, -2]
        ibase = 0

        for a_nn in a_lat:
            lat = fcc.LatFccBase4(l_lat, a_nn)
            lnn = ln.LatNeighbor(lat)

            mat_type = lnn.typ
            fcc_dist = a_nn*fcc_d

            # check that all types are of specie 0
            self.assertTrue(la.norm(mat_type - np.zeros_like(mat_type)) == 0)

            # check distances
            for ind in range(4):
                self.assertTrue(np.allclose(fcc_dist[ind], lnn.u_dist[ibase][ind+1]))

            # check multiplicity
            for ind in range(4):
                for ind_p in lst:
                    neighbor = lnn.neighbor(ind_p, ind+1)
                    self.assertTrue(len(neighbor) == fcc_mult[ind])

    def test_mgzn2(self):

        gam = np.sqrt(2 / 3)
        # the first distances of nearest neighbors are
        mgzn2_d = np.array([[np.sqrt(11/12), 1, 3/2], [gam, np.sqrt(11/12), np.sqrt(2)]])
        # with multiplicity
        mgzn2_mult = np.array([[12, 4, 15], [6, 6, 12]])

        l_lat = 3
        # list of lattice constants
        a_lat = np.array([0.5, 1.0, 3.0])
        # a_lat = np.array([1])
        # list of points
        lst = [0, 1, -1, -2]

        for a_nn in a_lat:
            lat = mgzn.LatMgZn2Base12(l_lat, a_nn, gam)
            lnn = ln.LatNeighbor(lat)

            mgzn2_dist = a_nn * mgzn2_d

            # check distances
            for ind in range(3):
                for ind_p in lst:
                    itype = lnn.typ[ind_p]
                    ibase = lnn.base[ind_p]
                    self.assertTrue(np.allclose(mgzn2_dist[itype][ind], lnn.u_dist[ibase][ind + 1]))

            # check multiplicity
            for ind in range(3):
                for ind_p in lst:
                    typ = lnn.typ[ind_p]
                    neighbor = lnn.neighbor(ind_p, ind + 1)
                    self.assertTrue(len(neighbor) == mgzn2_mult[typ][ind])


if __name__ == '__main__':
    unittest.main()
