"""
:module: UnitTest_ForcefieldReader.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the force field readers

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> April 2019
.. history:           
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - Added more tests for more general force field reader functionality
..                  - Added more tests for new forcefield reader functionality
..
..                Alex Upah <alexupah@iastate.edu> June 2022
..                  - Removed nonfunctional tests
..                  - Added tests for new forcefield reader functionality
"""

import unittest
import numpy as np
from hoodlt.Data.Forcefield import ForceFieldReader as ffrdr


class TestNameParser(unittest.TestCase):

    def setUp(self):
        self.ff_reader = ffrdr.ForceFieldReader('ncs-in-solvent')
        self.ff_reader_special = ffrdr.ForceFieldReader('ncs-in-solvent',)

    def test_get_molecular_weight(self):
        self.assertAlmostEqual(196.96657, self.ff_reader.get_molecular_weight('Au'))

    def test_get_charge(self):
        self.assertAlmostEqual(-.6705, self.ff_reader.get_charge('NH2'))
        self.assertEqual(self.ff_reader.get_charge('CH3'), 0.0)

    def test_get_nbnd_sigma_single_particle(self):
        self.assertEqual(1.559, self.ff_reader.get_nbnd_sigma_single_particle('HA'))

    def test_get_list(self):
        self.assertEqual(self.ff_reader.get_list_attributes(), ['groups', 'nonbonded', 'bond', 'angle', 'dihedral'])
        self.assertEqual(self.ff_reader.get_list_attributes_without_cutoff(), ['bond', 'angle', 'dihedral'])
        
    def test_has_rcut(self):
        self.assertEqual(self.ff_reader.has_rcut('bond', 'harmonic'), False)
        self.assertEqual(self.ff_reader.has_rcut('angle', 'harmonic'), False)
        self.assertEqual(self.ff_reader.has_rcut('dihedral', 'opls'), False)
        
    def test_get_list_attributes(self):
        self.assertEqual(self.ff_reader.get_list_attributes(), ['groups', 'nonbonded', 'bond', 'angle', 'dihedral'])
        
    def test_get_non_bonded(self):
        self.assertEqual(self.ff_reader.get_non_bonded('lj', 'CH2', 'CH3'),
                         {'sigma': 3.905, 'epsilon': 0.006233814879364159})
        self.assertEqual(self.ff_reader.get_non_bonded('lj', 'CH2', 'S'),
                         {'sigma': 4.1775, 'epsilon': 0.0074550613129337785})
    
    def test_get_bond_r0(self):
        self.assertEqual(self.ff_reader.get_bond_r0('S-CH2'), 1.82)
        
    def test_get_angle_t0(self):
        self.assertEqual(self.ff_reader.get_angle_t0('S-CH2-CH2', pot='harmonic'), 1.996)
        
    def test_get_potentials(self):
        self.assertEqual(self.ff_reader.get_potentials('nonbonded'), ['lj'])
        self.assertEqual(self.ff_reader.get_potentials('bond'), ['harmonic'])
        self.assertEqual(self.ff_reader.get_potentials('angle'), ['harmonic'])
        self.assertEqual(self.ff_reader.get_potentials('dihedral'), ['opls'])
    
    def test_get_potentials_params(self):
        val0 = self.ff_reader.get_potentials_params('nonbonded', 'lj', 'Au')
        val1 = self.ff_reader.get_potentials_params('bond', 'harmonic', 'CH2-CH2')
        val2 = self.ff_reader.get_potentials_params('angle', 'harmonic', 'CH2-CH2-CH2')
        val3 = self.ff_reader.get_potentials_params('dihedral', 'opls', 'CH2-CH2-CH2-CH2')
        
        self.assertEqual(val0, {'sigma': 1, 'epsilon': 2.3266791})
        self.assertEqual(val1, {'k': 25.5076414932, 'r0': 1.53})
        self.assertEqual(val2, {'k': 5.3859174233, 't0': 1.955})
        self.assertEqual(val3, {'k1': 0.21280861669588, 'k2': -0.127006165471831,
                                'k3': 0.228923023648524, 'k4': -0.0589895798559402})
      
    def test_get_potentials_dict(self):
        self.assertEqual(self.ff_reader.get_potentials_dict('nonbonded', 'lj'),
                         {'sigma': 'length', 'epsilon': 'energy'})
        self.assertEqual(self.ff_reader.get_potentials_dict('bond', 'harmonic'), {'k': 'bond', 'r0': 'length'})
        
    def test_get_mixing_rules(self):
        self.assertEqual(self.ff_reader._get_mixing_rules('lj'), {'sigma': 'arithmetic', 'epsilon': 'geometric'})


if __name__ == '__main__':
    unittest.main()
