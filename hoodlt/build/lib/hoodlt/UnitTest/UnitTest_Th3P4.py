"""
:module: UnitTest_Th3P4.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the Th3P4 unit class

.. moduleauthor:: Alex Travesset, December 2019
"""


import numpy as np
import numpy.linalg as la
import unittest
import hoodlt.Lattices.Lat2.Th3P4_lat as th3p4


class TestGeneralLJ(unittest.TestCase):

    def test_pf(self):
        """Test the pf
        """

        print('Test1: Packing fraction in Th3P4')

        gam = np.array([0.2, 0.288888888889, 0.377777777778, 0.466666666667, 0.555555555556, 0.644444444444])
        pf = np.array([0.649695566463, 0.663503516654, 0.6890501057, 0.729947222478, 0.680255514375, 0.635896570495])

        num = gam.shape[0]

        for ind in range(num):
            lat = th3p4.LatTh3P4Base28(3, 1, gam[ind])
            pf_c = lat.pf()
            self.assertTrue(np.allclose(pf_c, pf[ind], atol=0.001))

    def test_pnts(self):
        """"Test that the primitive vectors are consistent with their sizes
        """

        eps = 1e-5

        print('Test2: Points are OK in Th3P4')

        gam = np.array([0.4, 0.45, 0.55, 0.575, 0.6, 0.65, 0.665, 0.7, 0.775, 0.8])

        num = gam.shape[0]

        for ind in range(num):
            lat = th3p4.LatTh3P4Base28(3, 1, gam[ind])
            dd = np.array([1.0, gam[ind]])
            for ind1 in range(2):
                for ind2 in range(2):
                    ma = range(12 + ind1*4)
                    mb = range(12 + ind2*4)
                    dcl = 0.5 * (dd[ind1] + dd[ind2])
                    for ind3 in ma:
                        for ind4 in mb:
                            v1 = lat.get_v([ind1, ind3])
                            v2 = lat.get_v([ind2, ind4])
                            dst = la.norm(v1-v2)
                            if dst > eps:
                                self.assertTrue((dst > dcl - eps))

    def test_nei(self):

        print('checking neighbors')
        eps = 1e-7
        gam = 0.6
        lat = th3p4.LatTh3P4Base28(3, 1, gam)

        v0 = lat.get_v([0, 0])

        num_neighs = 8
        nn = np.array([12, 19, 22, 13, 23, 25, 26, 16])

        va = np.zeros([num_neighs, 3])
        for ind in range(num_neighs):
            inp = nn[ind]
            va[ind] = lat.get_v([0, inp])

        va[0] = va[0] - lat.get_a(1) - lat.get_a(2)
        va[1] = va[1] - lat.get_a(1)
        va[2] = va[2] - lat.get_a(1)
        va[6] = va[6] + lat.get_a(0) - lat.get_a(1)
        va[7] = va[7] + lat.get_a(0)

        dd = 0.5*(1+gam)
        for ind in range(num_neighs):
            dst = la.norm(v0-va[ind])
            self.assertTrue(np.abs(dd-dst) < eps)

        gam = 0.4

        lat = th3p4.LatTh3P4Base28(3, 1, gam)
        v0 = lat.get_v([0, 0])
        num_neighs = 8
        nn = np.array([3, 10, 8, 11, 2, 5, 4, 9])

        for ind in range(num_neighs):
            inp = nn[ind]
            va[ind] = lat.get_v([0, inp])

        va[0] = va[0] - lat.get_a(1)
        va[1] = va[1] - lat.get_a(1)
        va[4] = va[4] + lat.get_a(0) - lat.get_a(1)
        va[5] = va[5] + lat.get_a(0) - lat.get_a(1)
        va[6] = va[6] + lat.get_a(0) - lat.get_a(2)
        va[7] = va[7] + lat.get_a(0)

        for ind in range(num_neighs):
            dst = la.norm(v0-va[ind])
            self.assertTrue(np.abs(1-dst) < eps)


if __name__ == '__main__':
   unittest.main()
