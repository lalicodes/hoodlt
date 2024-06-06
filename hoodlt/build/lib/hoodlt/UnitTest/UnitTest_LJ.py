"""
:module: UnitTest_LJ

:platform: Unix, Windows
:synopsis: Defines the unit test for the Lennard Jones model

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2015
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test format
"""


import numpy as np
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as Fcc
import hoodlt.Lattices.Lat1.hcp_lat_Mixt as Hcp
import hoodlt.Potentials.overrp_pot_cutoff as ov
import hoodlt.Potentials.general_LJ as Lj
import hoodlt.D_matrix.ThermFunctions as ThFs
import unittest


class TestLJ(unittest.TestCase):

    def setUp(self):
        num_points = 12
        alat = 1
        cut_off = 4.1
        p_exp = 12
        q_exp = 6
        fcc_b = Fcc.LatFccBase4(num_points, alat)
            
        self.pot12 = ov.OverrpCut(p_exp, 1.0, 1.0, cut_off)
        self.d12 = ThFs.ThermoFunction(fcc_b, self.pot12)
            
        self.pot6 = ov.OverrpCut(q_exp, 1.0, 1.0, cut_off)
        self.d6 = ThFs.ThermoFunction(fcc_b, self.pot6)
            
        self.vlj = Lj.LJGen(p_exp, q_exp, 1.0, 1.0, cut_off)
        self.d_lj = ThFs.ThermoFunction(fcc_b, self.vlj)
        
    def test_energy(self):
        e_lj = self.d_lj.energy()
        e12 = self.d12.energy()
        e6 = self.d6.energy()
        elj_calc = 4 * (e12 - e6)
        self.assertAlmostEqual(e_lj, elj_calc, 12)
        
    def test_enthalpic_pressure(self):
        p_lj = self.d_lj.prs_enthalpic()
        p12 = self.d12.prs_enthalpic()
        p6 = self.d6.prs_enthalpic()
        plj_calc = 4 * (p12 - p6)

        self.assertAlmostEqual(p_lj, plj_calc, 8)
        
    def test_dmatrix(self):
        val = np.sum(np.abs(self.d_lj.mat_d - 4*(self.d12.mat_d - self.d6.mat_d)))
        self.assertAlmostEqual(val, 0.0)
        
        
if __name__ == '__main__':
    unittest.main()   
