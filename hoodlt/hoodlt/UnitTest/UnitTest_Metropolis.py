"""
:module: UnitTest_Metropolis.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the Metropolis class

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> , April 2017
"""

import hoodlt.MonteCarlo.Metropolis as Metro
import numpy as np
import unittest


class TestMetropolis(unittest.TestCase):

    def test_Pass(self):
        passed = 0
        times = 1000
        for i in range(0, times):
            delta_e = np.random.rand()
            rnd = np.exp(-(delta_e + 10e-8))
            if Metro.metropolis_test(delta_e, rnd):
                passed = passed+1
        self.assertEqual(passed, times)

    def test_Fail(self):
        passed = 0
        times = 1000
        for i in range(0, times):
            delta_e = np.random.rand()
            rnd = np.exp(-(delta_e - 10e-8))
            if Metro.metropolis_test(delta_e, rnd):
                passed = passed+1
        self.assertEqual(passed, 0)

    def test_Fail_Same(self):
        passed = 0
        times = 1000
        for i in range(0, times):
            delta_e = np.random.rand()
            rnd = np.exp(-delta_e)
            if Metro.metropolis_test(delta_e, rnd):
                passed += 1
        self.assertEqual(passed, 0)

if __name__ == '__main__':
    unittest.main()
