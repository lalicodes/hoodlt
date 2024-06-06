"""
:module: UnitTest_ptp_335.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the icosahedral group class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> , May 2017
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -removed unneccessary print statements
"""

import numpy as np
import unittest
import hoodlt.Groups.Quaternion as Qt
import hoodlt.Groups.Polytope_335 as PtP


class TestPtP(unittest.TestCase):

    def test_group_sum(self):
        """tests sums
        """

        yg = PtP.Polytope335()
        for ind1 in np.random.randint(0, yg.ico.order, size=20):
            qico = yg.ico.element_indx(ind1)
            for ind2 in np.random.randint(0, yg.num_vertices, size=8):
                qref = Qt.Quaternion(yg.vertices.quat[ind2])
                element_in_group = False
                for ind3 in range(yg.num_vertices):
                    qg = Qt.Quaternion(yg.vertices.quat[ind3])
                    if qg.approx(qref*qico):
                        element_in_group = True
                self.assertTrue(element_in_group)

if __name__ == '__main__':
    unittest.main()
