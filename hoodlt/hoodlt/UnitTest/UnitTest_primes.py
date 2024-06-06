"""
:module: UnitTest_primes.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the calculation of primes

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February 2016
"""

import hoodlt.Utils.Prime_Factors as Pf
import unittest


class TestPrimes(unittest.TestCase):

    def test_primes8(self):
        """Test list of primes for 8
        """
        u_8 = list(Pf.divisor_generator(8))
        u_8_exact = [1, 2, 4, 8]

        self.assertTrue(u_8, u_8_exact)

    def test_primes50(self):
        """Test list of primes for 50
        """
        u_50 = list(Pf.divisor_generator(50))
        u_50_exact = [1, 2, 5, 10, 25, 50]

        self.assertTrue(u_50, u_50_exact)

if __name__ == '__main__':
    unittest.main()
