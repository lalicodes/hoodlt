"""
:module: Dmatrix_Mixt_potential_decorator
:platform: Unix, Windows
:synopsis: Defines the decorator to call the potential function

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December2015
"""

import numpy as np


def compute_norm(func):
    """Defines the decorator for the argument function. This decorator does not keep track of the original matrix

    :param func: function to be wrapped
    :return: wrapper
    :rtype: method
    """
    def wrapper(self, r):
        """Defines the wrapper

        :param: r vector distance
        :return: wrapper function
        :rtype: method
        """

        rr = np.sum(r*r, axis=0)

        is_an_array = isinstance(rr, np.ndarray)

        if is_an_array:
            rrmat = rr[np.logical_and(rr > self.epsilon*self.epsilon, rr < self.cut*self.cut)]
            if rrmat.size == 0:
                return np.array((0.0,))
        else:
            if rr < self.cut * self.cut:
                rrmat = rr
            else:
                return 0.0

        return func(self, np.sqrt(rrmat))

    return wrapper


def compute_norm_keep_matrix(func):
    """Returns the value of the function as a numpy array with the same structure as the input array.
    It is slower than compute_norm but needs to be considered in the calculation of the D-matrix

    :param: func
    :return: wrapper
    :rtype: method
    """

    def wrapper(self, *args):
        """Defines the wrapper

        :param: r vector distance
        :return: wrapper function
        :rtype: method
        """

        r = args[0]
        rr = np.sum(r*r, axis=0)

        is_an_array = isinstance(rr, np.ndarray)

        if is_an_array:
            rr_ind = np.where(np.logical_and(rr > self.epsilon*self.epsilon, rr < self.cut*self.cut))
            rrmat = rr[rr_ind]
            if rrmat.size == 0:
                return np.array((0.0,))
        else:
            if rr < self.cut * self.cut:
                rrmat = rr
            else:
                return 0.0

        val = func(self, np.sqrt(rrmat))

        if is_an_array:
            res = np.zeros_like(rr)
            res[rr_ind] = val
            val = res

        return val

    return wrapper
