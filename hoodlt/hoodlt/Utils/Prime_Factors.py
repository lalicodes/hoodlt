"""
:module: Prime_Factors
:platform: Unix, Windows
:synopsis: Utility to calculate prime factors

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2016
"""

import numpy as np


def divisor_generator(n):
    """Returns the prime factors up to value n
    
    :param n: maximum value
    :return: all prime factors
    :rtype: generator
    """
    large_divisors = []
    for i in range(1, int(np.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor
