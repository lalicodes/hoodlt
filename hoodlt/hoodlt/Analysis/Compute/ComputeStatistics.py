"""
:module: ComputeStatistics
:Platform: Windows, Unix
:synopsis: Statistical tools for analyzing time series

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May 2022
.. history:
..
"""

import numpy as np
from scipy.stats import sem


class StatComputeTimeSeries(object):
    """
    class to analyze time series

    """

    def __init__(self, mat, ini_val=0, end_val=None):
        """Get the standard deviation of a given quantity

            :param mat: numpy array
            :param ini_val: initial index to consider
            :param end_val: final index to consider
        """

        self.mat = mat

        if end_val is None:
            self.end_val = mat.shape[0]
        else:
            self.end_val = end_val

        self.ini_val = ini_val

        self.num_points = self.end_val - self.ini_val

    def compute_average(self):
        """Get the standard deviation of a given quantity

            :return: ndarray
        """

        return np.average(self.mat[self.ini_val: self.end_val], axis=0)

    def compute_std_error(self):
        """Get the standard deviation of a given quantity

            :return: ndarray
        """

        return sem(self.mat[self.ini_val: self.end_val], axis=0)

    def compute_std_deviation(self):
        """Get the standard deviation of a given quantity

            :return: ndarray
        """

        return np.std(self.mat[self.ini_val: self.end_val], axis=0)

    def compute_autocorrelation(self, max_num=20):
        """Get the correlation in time series

            :param max_num: maximum time step to compute correlation
            :return: ndarray
        """

        mt = self.mat[self.ini_val:self.end_val]

        r1 = [1.0 if ind == 0 else np.corrcoef(mt[ind:], mt[:-ind])[0, 1] for ind in range(max_num)]

        return np.array(r1)

    def compute_std_error_by_block(self):
        """Estimation of the time correlation

            :return: ndarray
        """

        k_val = int(np.floor(np.log(self.num_points)/np.log(2.0)))

        new_ini_val = self.ini_val+self.num_points-2**k_val
        mt = self.mat[new_ini_val:self.end_val]
        var_egy = np.zeros(k_val - 1)
        for indl in range(k_val-1):
            mat_blk = np.zeros(2**(k_val-indl))
            for indp in range(2 ** (k_val - indl)):
                i_val = indp * 2 ** indl
                e_val = (indp + 1) * 2 ** indl
                mat_blk[indp] = np.sum(mt[i_val:e_val]) / 2 ** indl
            var_egy[indl] = np.var(mat_blk, ddof=1) / 2 ** (k_val - indl)

        return np.sqrt(var_egy)

    def compute_thermalization(self):
        """
        fit to linear data to establish if the data has thermalized

        :return
        """
        x = np.arange(self.ini_val, self.end_val)
        y = self.mat[self.ini_val:self.end_val]

        return x, np.polyfit(x, y, 1)
