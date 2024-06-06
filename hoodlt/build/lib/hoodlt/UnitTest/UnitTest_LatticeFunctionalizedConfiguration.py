"""
:module: UnitTest_LatticeFunctionalizedConfiguration.py
:platform: Unix, Windows
:synopsis: Defines the unit test for LatticeFunctionalizedConfiguration class

.. moduleauthor:: Xun Zha <xzha@iastate.edu> March 2021
.. history:
..                Alex Upah <alexupah@iastate.edu> June 2022
..                  - Updated tests to properly function following hooldt updates
"""


import unittest
import numpy as np
import hoodlt.Data.Modelconfigurations.LatticeFunctionalizedConfiguration as LatConfig
from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
import hoodlt.Lattices.Lat1.fcc_lat_Mixt as Fcc
import hoodlt.Lattices.Lat2.MgZn2_lat as MgZn2
from hoodlt.Data.Modelnanoparticles.AuS import AuS
from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand


class TestLatFCC(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # forcefield name
        forcefield = "dry-ncs"

        lat = Fcc.LatFccBase4(l_value=2, a_nn_e=50)
        builder = ConfigurationBuilder()
        core = AuS(forcefield, 201, 1)
        lig = [HydrocarbonLigand(repeats=11, ff=forcefield)] * core.graft_num
        builder.build_lattice([core], [lig], lat)
        cls.conf = builder.conf

    def test_box(self):
        self.assertTrue(np.allclose(self.conf.box_l, np.concatenate((np.ones(3)*100*np.sqrt(2), np.zeros(3)))))

    def test_rotation(self):
        self.assertTrue(np.allclose(self.conf.particles[0].core.orientation[0], np.array([1, 0, 0, 0])))

    def test_bond(self):
        bond_names = ['CTR-CTR1', 'CTR-CTR2']
        
        self.conf.add_bonds({'CTR-CTR1': (0, 1), 'CTR-CTR2': (0, 2)})
        self.assertListEqual(self.conf.bonds_types, bond_names)
        self.assertEqual(len(self.conf.bonds), 2)
        self.assertEqual(len(self.conf.bonds[0]), 192)
        self.assertEqual(len(self.conf.bonds[1]), 48)


class TestLatMgZn2(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # forcefield name
        forcefield = "dry-ncs"

        lat = MgZn2.LatMgZn2Base12(l_value=2, a_nn_e=50, gamma=1)
        builder = ConfigurationBuilder()
        core_a = AuS(forcefield, 1072, 2)
        core_b = AuS(forcefield, 201, 1)
        lig_a = [HydrocarbonLigand(repeats=11, ff=forcefield)] * core_a.graft_num
        lig_b = [HydrocarbonLigand(repeats=11, ff=forcefield)] * core_b.graft_num
        builder.build_lattice([core_a, core_b], [lig_a, lig_b], lat)
        cls.conf = builder.conf

    def test_box(self):
        self.assertTrue(np.allclose(self.conf.box_l, np.concatenate((np.array([2, np.sqrt(3), np.sqrt(8/3.)])*100,
                                                                     np.array([-1/np.sqrt(3), 0, 0])))))

    def test_rotation(self):
        self.assertTrue(np.allclose(self.conf.particles[0].core.orientation[0], np.array([1, 0, 0, 0])))

    def test_bond(self):
        bond_names = ['CTR-CTR1', 'CTR-CTR2', 'CTR-CTR3']
        self.conf.add_bonds({'CTR-CTR1': (0, 1), 'CTR-CTR2': (0, 2), 'CTR-CTR3': (1, 1)})
        self.assertListEqual(self.conf.bonds_types, bond_names)
        self.assertEqual(len(self.conf.bonds), 3)
        self.assertEqual(len(self.conf.bonds[0]), 192)
        self.assertEqual(len(self.conf.bonds[1]), 32)
        self.assertEqual(len(self.conf.bonds[2]), 96)


class TestBinaryLatticeAddingBonds(unittest.TestCase):

    def setUp(self):
        # forcefield name
        forcefield = "dry-ncs"

        lat = MgZn2.LatMgZn2Base12(l_value=2, a_nn_e=50, gamma=1)
        builder = ConfigurationBuilder()
        core_a = AuS(forcefield, 1072, 2)
        core_b = AuS(forcefield, 201, 1)
        lig_a = [HydrocarbonLigand(repeats=11, ff=forcefield)]*core_a.graft_num
        lig_b = [HydrocarbonLigand(repeats=11, ff=forcefield)]*core_b.graft_num
        builder.build_lattice([core_a, core_b], [lig_a, lig_b], lat)
        self.conf = builder.conf

    def test_bond(self):
        self.conf.add_bonds({'CTR-CTR1': (1, 1)})
        self.assertListEqual(self.conf.bonds_types, ['CTR-CTR1'])
        self.assertEqual(len(self.conf.bonds), 1)
        self.assertEqual(len(self.conf.bonds[0]), 96)
        self.assertListEqual(self.conf.bonds, [self.conf.bonds[0]])

        self.conf.add_bonds({'CTR-CTR2': (1, 1)})
        self.assertListEqual(self.conf.bonds_types, ['CTR-CTR2'])
        self.assertEqual(len(self.conf.bonds), 1)
        self.assertEqual(len(self.conf.bonds[0]), 96)

    def test_error1(self):
        try:
            self.conf.add_bonds({'CTR-CTR1': (0, 5.378)})
        except IndexError:
            pass
        else:
            raise AssertionError('IndexError was not raised.')

    def test_error2(self):
        try:
            self.conf.add_bonds({'CTR-CTR1': (1, -7)})
        except IndexError:
            pass
        else:
            raise AssertionError('IndexError was not raised.')


class TestLatRotated(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # forcefield name
        forcefield = "dry-ncs"

        lat = Fcc.LatFccBase4(l_value=2, a_nn_e=50)
        builder = ConfigurationBuilder()
        core = AuS(forcefield, 201, 1)
        lig = [HydrocarbonLigand(repeats=11, ff=forcefield)]*core.graft_num
        ornt = np.zeros((32, 4))
        ornt[:] = np.array([1, 0, 0, 0])
        builder.build_lattice([core], [lig], lat)
        cls.conf = builder.conf

    def test_rotation(self):
        self.assertTrue(np.allclose(self.conf.particles[0].core.orientation[0], np.array([1, 0, 0, 0])))


if __name__ == '__main__':
    unittest.main()
