"""
:module: UnitTest_Redefine_lattice.py
:platform: Unix, Windows
:synopsis: testing the ability to redefine lattices by changing the lattice constant, which is a new addition.

.. moduleauthor:: Alex Travesset, February 2020
"""


import numpy as np
import unittest
import hoodlt.Lattices.Lat2.CaTiO3b_lat as cat3b


class TestGeneral(unittest.TestCase):

    def test_pf(self):
        """Test the pf
        """

        a_nn = 1.0
        gam = 0.4

        lat1 = cat3b.LatCaTiO3Base5(3, a_nn, gam)
        lat2 = cat3b.LatCaTiO3Base5(3, a_nn, gam)

        a_nn_new = 3.0
        fac = a_nn_new/a_nn

        lat2.redefine_lattice_constant(a_nn_new)

        # test primitive vectors
        for ind in range(3):
            self.assertTrue(np.allclose(fac*lat1.get_a(ind), lat2.get_a(ind)))

        # test basis
        for ind1 in range(2):
            for ind2 in range(lat1.typ[ind1]):
                self.assertTrue(np.allclose(fac*lat1.get_v([ind1, ind2]), lat2.get_v([ind1, ind2])))

        # test box
        self.assertTrue(np.allclose(lat1.l_box()*fac, lat2.l_box()))

        # test volume
        self.assertTrue(np.allclose(lat1.vol_unit_cell()*fac**3, lat2.vol_unit_cell()))

        # test packing fraction
        self.assertTrue(np.allclose(lat1.pf()/fac**3, lat2.pf()))

        # test number of points
        self.assertTrue(np.allclose(lat1.num_pnts(), lat2.num_pnts()))


if __name__ == '__main__':
    unittest.main()
