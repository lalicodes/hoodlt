"""
:module: UnitTest_HydrocarbonLigand.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the HydrocarbonLigand class

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> , April 2017
.. history:
..                Alex Upah <alexupah@iastate.edu> June 2022
..                   -removed nonfunctional tests
"""

from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand
import unittest
import numpy as np
import os
import shutil
import copy as cp


class TestHydrocarbon(unittest.TestCase):

    def setUp(self):
        self.dir = "dumpfolder/"
        directory = os.path.dirname(self.dir)

        self.ff = "opls_dry-ncs_mix"

        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_Positions(self):
        chain = HydrocarbonLigand(11, self.ff)
        self.assertTrue(np.allclose(chain.position[0], np.array([0.0, 0.0, 0.0])))
        self.assertAlmostEqual(chain.position[12][0], 15.273590175884479)

    def test_shift(self):
        chain = HydrocarbonLigand(2, self.ff)
        chain2 = HydrocarbonLigand(2, self.ff)
        chain.shift([1, 1, 1])
        for i in range(0, 4):
            self.assertTrue(np.allclose(chain.position[i], np.add([1, 1, 1], chain2.position[i])))
        self.assertTrue(np.allclose(chain.position[0], [1, 1, 1]))

    def test_bonds(self):
        chain = HydrocarbonLigand(8, self.ff)
        # print chain.get_bonds()
        self.assertEqual(len(chain.get_bonds()), 9)

    def test_angles(self):
        chain = HydrocarbonLigand(8, self.ff)
        # print chain.get_angles()
        self.assertEqual(len(chain.get_angles()), 8)

    def test_dihedrals(self):
        chain = HydrocarbonLigand(8, self.ff)
        # print chain.get_dihedrals()
        self.assertEqual(len(chain.get_dihedrals()), 7)

    def test_bond_types(self):
        chain = HydrocarbonLigand(11, self.ff)
        # print chain.get_dihedrals()
        self.assertEqual(chain.get_bond_types()[0], 'S-CH2')
        self.assertEqual(len(chain.get_bond_types()), 12)

    def test_angle_types(self):
        chain = HydrocarbonLigand(11, self.ff)
        # print chain.get_dihedrals()
        self.assertEqual(chain.get_angle_types()[0], 'S-CH2-CH2')
        self.assertEqual(len(chain.get_angle_types()), 11)

    def test_dihedral_types(self):
        chain = HydrocarbonLigand(11, self.ff)
        self.assertEqual(chain.get_dihedral_types()[0], 'CH2-CH2-CH2-CH2')
        self.assertEqual(len(chain.get_dihedral_types()), 10)

    def test_mass(self):
        chain = HydrocarbonLigand(11, self.ff)
        self.assertEqual(33.0729, chain.mass[0])
        self.assertEqual(15.0345, chain.mass[12])
        for i in range(1, 12):
            self.assertEqual(14.0266, chain.mass[i])

    def test_orient_consistency(self):
        lig = HydrocarbonLigand(11, self.ff)
        part = cp.deepcopy(lig)
        lig.rotate([0, 0, 0], [0.5, 0.5, 0.5, 0.5])\
        
        self.assertFalse(np.allclose(lig.position, part.position))
        
        lig.rotate([0, 0, 0], [0.5, -0.5, -0.5, -0.5])
        self.assertTrue(np.allclose(lig.position, part.position))
        
    def test_scale(self):
        lig = HydrocarbonLigand(11, self.ff)
        part = cp.deepcopy(lig)

        lig.scale(4)

        self.assertTrue(np.allclose(lig.position, part.position * 4))


if __name__ == '__main__':
    unittest.main()
