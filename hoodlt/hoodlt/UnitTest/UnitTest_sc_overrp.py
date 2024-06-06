"""
:module: UNitTest_sc_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the classes implementing the sc lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.sc_lat_Mixt as Sc
import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import hoodlt.Potentials.overrp_pot_cutoff_Mixt as Ov
import unittest


class Testsc(unittest.TestCase):
    
    def setUp(self):
        self.sc = Sc.LatSC(20, 1.0)
        
        # these values are from Ashcroft and Mermin
        self.mat_p = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,  16])
        self.mat_sc = np.array([16.53, 10.38, 8.40, 7.47, 6.95, 6.63, 6.43, 6.29, 6.20, 6.14, 6.10, 6.07, 6.05])
        self.mat_err = np.array([1.4, 0.09, 0.01, 0.005, 0.005, 0.004, 0.004, 0.004, 0.004, 0.002, 0.002, 0.002, 0.002])
        self.cut_off = 9.1
        
    def test_lattice(self):
        for i in range(self.mat_p.size):
            vp = Ov.OverrpCutMixt(self.mat_p[i], np.array([1]), np.array([[1]]), self.cut_off)
            val = 2 * Dm.DMatrix(self.sc, vp).energy()
            
            self.assertTrue(np.abs(self.mat_sc[i] - val) < self.mat_err[i])
        
    def test_pf(self):
        self.assertAlmostEqual(self.sc.pf(), np.pi/6.0, 12)
        

if __name__ == '__main__':
    unittest.main()
