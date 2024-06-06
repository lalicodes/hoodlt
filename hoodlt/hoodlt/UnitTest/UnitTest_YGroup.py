"""
:module: UnitTest_YGroup.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the icosahedral group class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> , May 2017
"""


import hoodlt.Groups.Yico_Group as YG
import numpy as np
import unittest


class TestIco(unittest.TestCase):

    def test_group_sum(self):
        """tests sums
        """

        yg = YG.YGroup()
        for ind1 in np.random.randint(0, yg.order, size=20):
            qval = yg.element_indx(ind1)
            for ind2 in np.random.randint(0, yg.order, size=8):
                qref = yg.element_indx(ind2)
                element_in_group = False
                for ind3 in range(yg.order):
                    qg = yg.element_indx(ind3)
                    if qg.approx(qval*qref):
                        element_in_group = True
                self.assertTrue(element_in_group)

if __name__ == '__main__':
    unittest.main()
