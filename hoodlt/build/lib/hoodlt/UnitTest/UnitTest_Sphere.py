"""
:module: UnitTest_Sphere.py
:platform: Unix, Windows
:synopsis: Defines the unit test for Sphere

.. moduleauthor:: Nathan Horst, April 2017
.. history:
..                Alex Upah <alexupah@iastate.edu> June 2022
..                  -removed non functional tests, updated function parameters
"""

import hoodlt.Data.Modelnanoparticles.Sphere as Sph
import numpy as np
import unittest
import copy as cp
import hoodlt.Data.Modelnanoparticles.SphereInfo as Si
import os
import shutil


class TestSphere(unittest.TestCase):

    def setUp(self):
        self.ff = "cg_atacticPS"
        self.dir = "dumpfolder/"
        directory = os.path.dirname(self.dir)
        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_orienter_consistency(self):

        sphere_info = Si.SphereInfo(10, 10, 100, 80)

        particle = Sph.Sphere(sphere_info, self.ff)
        part = cp.deepcopy(particle)
        
        particle.rotate([0.5, 0.5, 0.5, 0.5])
        particle.rotate([0.5, -0.5, -0.5, -0.5])
        self.assertTrue(np.allclose(particle.position, part.position))


if __name__ == '__main__':
    unittest.main()
