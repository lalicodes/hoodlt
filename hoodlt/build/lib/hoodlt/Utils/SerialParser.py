"""
:module: SerialParser
:platform: Unix, Windows
:synopsis: Parses .hist and .log files based on r0 values

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> May 2017
"""

from __future__ import division
import numpy as np


class SerialParser(object):
    """
        Parse .hist and .log files created by a serial run to be based on r0 value
    """
    def __init__(self, s_name, n_winds):
        """
        
        :param s_name: name of the serial run i.e. the file it was originally run from without the .gsd extension
        :param n_winds: the number of windows in the simulation
        """
        if type(s_name) == list:
            if len(s_name) > 1:
                raise ValueError("SerialParser only parses one file.")
            s_name = s_name[0]

        self.file = s_name
        self.n_winds = n_winds

    def make_histogram_data(self, len_unit=1):
        """
        
        :param len_unit: length unit to convert the data properly
        :return: histogram_data np.array((number of windows, points per window))
        """
        f = open(self.file + '.hist')
        datas = f.readlines()

        if len(datas) % self.n_winds != 0:
            raise ValueError('Illegal number of data points for ' + str(self.n_winds) + "r values.")

        points_per = int(len(datas) / self.n_winds)

        data_file = []

        for i in range(self.n_winds):
            data = [float(datas[ind]) for ind in range(points_per * i, points_per * (i + 1))]
            data_file.append(data)

        return np.array(data_file) * len_unit

    def make_log_data(self, column_index, e_unit=1):
        """
        :param column_index: index of the column in the log to make the data for
        :param conv_unit: multiplier for all the data
        :return:  np.array((number of windows, logs per window))
        """
        f = open(self.file + '.log')
        datas = f.readlines()

        datas = [data for data in datas if data.split()[0] != 'timestep']
        if len(datas) % self.n_winds != 0:
            raise ValueError('Illegal number of data points (' + str(len(datas)) + ') for ' + str(self.n_winds)
                             + "r values.")

        points_per = int(len(datas) / self.n_winds)

        data_file = []

        for i in range(self.n_winds):
            data = [float(datas[ind].split()[column_index]) for ind in range(points_per * i, points_per * (i + 1))]
            data_file.append(data)
        return e_unit * np.array(data_file)
