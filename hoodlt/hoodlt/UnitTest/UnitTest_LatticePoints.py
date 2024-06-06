"""
:module: UnitTest_LatticePoints.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the NanoInitGSD function

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> , February 2021
"""

import unittest
import numpy as np

import hoodlt.Utils.LatticePoints as Lp
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as fcc
import hoodlt.Lattices.Lat2.MgZn2_lat as mgzn


class TestLatticePoints(unittest.TestCase):

    def test_fcc(self):

        l_lat = 3
        # list of lattice constants
        a_lat = 1.0
        # list of points
        lat = fcc.LatFccBase4(l_lat, a_lat)

        lbl = Lp.LatPoints(lat)

        pnts = np.array([0, 5, 25, 50, 107], dtype='int')

        pnt_1 = np.array([0, 0, 0, 0, 0], dtype='int')
        pnt_2 = np.array([0, 1, 1, 0, 0], dtype='int')
        pnt_3 = np.array([0, 1, 0, 2, 0], dtype='int')
        pnt_4 = np.array([0, 2, 0, 1, 1], dtype='int')
        pnt_5 = np.array([0, 3, 2, 2, 2], dtype='int')

        pnts_sol = [pnt_1, pnt_2, pnt_3, pnt_4, pnt_5]

        # check that five coordinates are correct
        for ind, ind_n in enumerate(pnts):
            self.assertTrue(np.array_equal(lbl.coordinate_five(ind_n), pnts_sol[ind]))

        # check that single coordinates are correct
        for ind, ind_n in enumerate(pnts_sol):
            self.assertTrue(lbl.i_val(ind_n) == pnts[ind])

    def test_mgzn2(self):

        # gamma value
        gam = np.sqrt(2 / 3)
        # number of unit cells
        l_lat = 3
        # list of lattice constants
        a_lat = 1.0

        lat = mgzn.LatMgZn2Base12(l_lat, a_lat, gam)

        lbl = Lp.LatPoints(lat)

        pnts = np.array([0, 5, 25, 114, 201, 215], dtype='int')

        pnt_1 = np.array([0, 0, 0, 0, 0], dtype='int')
        pnt_2 = np.array([1, 1, 0, 0, 0], dtype='int')
        pnt_3 = np.array([0, 1, 2, 0, 0], dtype='int')
        pnt_4 = np.array([1, 2, 0, 0, 1], dtype='int')
        pnt_5 = np.array([1, 5, 1, 2, 1], dtype='int')
        pnt_6 = np.array([1, 7, 2, 2, 1], dtype='int')

        pnts_sol = [pnt_1, pnt_2, pnt_3, pnt_4, pnt_5, pnt_6]

        # check that five coordinates are correct
        for ind, ind_n in enumerate(pnts):
            self.assertTrue(np.array_equal(lbl.coordinate_five(ind_n), pnts_sol[ind]))

        # check that single coordinates are correct
        for ind, ind_n in enumerate(pnts_sol):
            self.assertTrue(lbl.i_val(ind_n) == pnts[ind])


if __name__ == '__main__':
    unittest.main()
