"""
:module: UnitTest_LJ_Mixt.py
:platform: Unix, Windows
:synopsis: Defines the unit test for general_LJ_Mixt

.. moduleauthor:: Alex Travesset, November 2018
"""


import numpy as np
import unittest
import os
import shutil
import hoodlt.Potentials.general_LJ_Mixt as lm
import hoodlt.Potentials.general_LJ as lj

class TestGeneralLJ(unittest.TestCase):

    def setUp(self):
        self.dir = "dumpfolder/"
        directory = os.path.dirname(self.dir)
        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_values(self):
        """Test values of the potential
        """

        p = 12
        q = 6
        sigma = np.array([1.0, 0.5])
        eps = np.ones([2, 2])
        cut_off = 4.0

        vals = [0.1, 1.0, 2.0, 1.5, np.pi]

        lj_mixt = lm.LJMixt(p, q, sigma, eps, cut_off)
        lj_mat = []

        for n in range(2):
            lj_mat.append([])
            for m in range(2):
                lj_mat[n].append(lj.LJGen(p, q, 0.5*(sigma[n]+sigma[m]), eps[n, m], cut_off))

        for x in vals:
            for n in range(2):
                for m in range(2):
                    self.assertTrue(np.allclose(lj_mat[n][m].der0(x), lj_mixt.pot[n][m].der0(x)))


if __name__ == '__main__':
    unittest.main()
