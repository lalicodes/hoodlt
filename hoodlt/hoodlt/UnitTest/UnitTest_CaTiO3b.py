"""
:module: UnitTest_CaTiO3b.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the CaTiO3b unit class

.. moduleauthor:: Alex Travesset, December 2019
"""


import numpy as np
import numpy.linalg as la
import unittest
import hoodlt.Lattices.Lat2.CaTiO3b_lat as cat3b


class TestGeneralLJ(unittest.TestCase):

    def test_pf(self):
        """Test the pf
        """

        print('Test1: Packing fraction in CaTiO3b')

        gam = np.array([0.2, 0.288888888889, 0.377777777778, 0.466666666667, 0.555555555556, 0.644444444444])
        pf = np.zeros_like(gam)

        def func(x):
            x1 = np.sqrt(2)-1
            x2 = 1/(2*np.sqrt(2)-1)
            if x < x1:
                return np.pi*(1+4*x**3)/6
            elif x < x2:
                return np.pi*np.sqrt(2)*(1+4*x**3)/(3*(1+x)**3)
            else:
                return np.pi*(1+4*x**3)/(48*x**3)

        for ind, val in enumerate(gam):
            pf[ind] = func(val)

        num = gam.shape[0]

        for ind in range(num):
            lat = cat3b.LatCaTiO3Base5(3, 1, gam[ind])
            pf_c = lat.pf()
            self.assertTrue(np.allclose(pf_c, pf[ind], atol=0.001))

    def test_pnts(self):
        """"Test that the primitive vectors are consistent with their sizes
        """

        eps = 1e-5

        gam = np.array([0.4, 0.45, 0.55, 0.575, 0.6, 0.65, 0.665, 0.7, 0.775, 0.8])

        num = gam.shape[0]

        for ind in range(num):
            lat = cat3b.LatCaTiO3Base5(3, 1, gam[ind])
            dd = np.array([1.0, gam[ind]])
            for ind1 in range(2):
                for ind2 in range(2):
                    ma = range(1 + ind1*3)
                    mb = range(1 + ind2*3)
                    dcl = 0.5 * (dd[ind1] + dd[ind2])
                    for ind3 in ma:
                        for ind4 in mb:
                            v1 = lat.get_v([ind1, ind3])
                            v2 = lat.get_v([ind2, ind4])
                            dst = la.norm(v1-v2)
                            if dst > eps:
                                self.assertTrue((dst > dcl - eps))


if __name__ == '__main__':
    unittest.main()
