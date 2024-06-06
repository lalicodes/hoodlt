"""
:module: UnitTest_SimParameters.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the SimParameters class

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu>, March 2017]
.. history:
..                Alex Upah <alexupah@iastate.edu>, June 2022
..                  -Updated to function with version three
"""


import hoodlt.HOOMD.SimParameters as Sparam
import numpy as np
import unittest


class TestSimParams(unittest.TestCase):

    def test_all(self):
        steps_sim = np.random.rand()
        steps_log = np.random.rand()
        steps_write = np.random.rand()
        steps_hist = np.random.rand()
        num_procs = np.random.rand()
        quantities_log = ['a', 'd', 'g']

        sp = Sparam.SimParameters(steps_sim, steps_log, quantities_log, steps_write, steps_hist, num_procs)
        self.assertEqual(steps_sim, sp.steps_sim)
        self.assertEqual(steps_log, sp.steps_log)
        self.assertEqual(steps_write, sp.steps_write)
        self.assertEqual(num_procs, sp.num_procs)
        self.assertEqual(steps_hist, sp.steps_hist)
        self.assertEqual(quantities_log[0], sp.quantities_log[0])
        self.assertEqual(quantities_log[1], sp.quantities_log[1])
        self.assertEqual(quantities_log[2], sp.quantities_log[2])


if __name__ == '__main__':
    unittest.main()
