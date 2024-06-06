"""
:module: UnitTest_Dmatrix_Pressure

:platform: Unix, Windows
:synopsis: Defines the unit test for the enthalpic pressure

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test to fit proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as Fcc
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import unittest


class TestDmatrixPressure(unittest.TestCase):
    
    def setUp(self):
        self.alat = [1.0, 0.5]
        self.cut_off = 2.1
        self.lsize = 6
        self.pval = 6

    def test_enthalpic(self):
    
        for i in range(len(self.alat)):
            lat = Fcc.LatFccBase4(self.lsize, self.alat[i])
            pot = Ov.OverrpCutMixt(self.pval, np.array([1.0]), np.array([[1.0]]), self.cut_off)
            dmat = Dm.DMatrix(lat, pot)
            egy = dmat.energy()
            prs_v = self.pval * egy / 3.0

            prs = dmat.prs_enthalpic()
            
            self.assertAlmostEqual(prs_v, prs, 12)


if __name__ == '__main__':
    unittest.main()
