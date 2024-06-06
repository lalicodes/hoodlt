"""
:module: SwapParser
:platform: Unix, Windows
:synopsis: Parses .hist and .log files based on r0 values

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> May 2017
"""

from __future__ import division
import numpy as np
import gsd.hoomd


class SwapParser(object):
    """
        Parse .hist and .log files created based on processor to be based on r0 value
    """
    def __init__(self, swap_files, swap=True):
        """
        
        :param swap_files: list of files containing the swapping information
        :param swap: whether or not any swaps were performed, defaults to true
        """

        if type(swap_files) == str:
            swap_files = [swap_files]

        self.num_procs = len(swap_files)
        self.info = []
        self.files = swap_files
        if swap:
            for ind, swap_file in enumerate(swap_files):
                f = open(swap_file + '.switch')
                data = f.readlines()
                if ind != 0:
                    temp = self.num_swaps
                else:
                    temp = (len(data) - 1)
                self.num_swaps = (len(data) - 1)
                if temp != self.num_swaps:
                    raise ValueError("All files must have the same number of swaps")
                self.info.append([[swap_file, tag[:-1]] for ind, tag in enumerate(data[:-1])])
        else:
            self.info = [[[the_file, the_file]] for the_file in swap_files]
            self.num_swaps = 1

    def make_histogram_data(self, len_unit=1):
        """
        :param len_unit: multiplier for all the distances/energies/whatever is being histogrammed
        :return: numpy array with the data
        """

        check = open(self.info[0][0][0] + '.hist', 'r')
        lines = len(check.readlines())
        data_file = [[] for i in range(self.num_procs)]
        for proc in self.info:
            f = open(proc[0][0] + '.hist', 'r')
            x = f.readlines()
            if len(x) != lines:
                raise ValueError("All files must be the same length. " + str(len(x)) + " doesn't equal " + str(lines))
            points_per = int(len(x) / self.num_swaps)
            for ind, y in enumerate(proc):
                for z in range(points_per):
                    data_file[self.files.index(y[1])].append(float(x[z + points_per * ind]))

        return len_unit * np.array(data_file)

    def make_log_data(self, column_index, e_unit=1):
        """
        
        :param column_index: log column to make the data for; indexing starts at 0
        :param conv_unit: multiplier for all the distances/energies/whatever is being logged
        :return: numpy array with the data
        """
        check = open(self.info[0][0][0] + '.log')
        x = check.readlines()
        points = len([data for data in x if data.split()[0] != 'timestep'])
        points_per = int(points / self.num_swaps)
        data_file = np.zeros((self.num_procs, points))

        for proc in self.info:
            f = open(proc[0][0] + '.log', 'r')
            x = f.readlines()
            x = [data for data in x if data.split()[0] != 'timestep']
            if len(x) % self.num_swaps != 0:
                raise ValueError("Illegal number of data points " + str(len(x)) + " is not divisible by " +
                                 str(self.num_swaps))
            for ind, y in enumerate(proc):
                for z in range(points_per):
                    data_point = x[z + points_per * ind].split()[column_index]
                    data_file[self.files.index(y[1])][z + points_per * ind] = float(data_point)

        return e_unit * data_file

    def parse_trajectories(self):
        """
        
        :return: 2d list of snapshots (numer of procs, number of frames per file) 
        """
        sys = gsd.hoomd.open(self.info[0][0][0] + '_dump.gsd', mode='rb')
        frames = len(sys)
        if frames % self.num_swaps != 0:
            raise ValueError("Illegal combination of frames, " + str(frames) + " and swaps, " + str(self.num_swaps))
        frames_per = int(frames / self.num_swaps)
        data_file = [[] for i in range(self.num_procs)]

        for proc in self.info:
            f = gsd.hoomd.open(proc[0][0] + '_dump.gsd', mode='rb')
            if len(f) != frames:
                raise ValueError("All files must have the same number of frames")
            for ind, y in enumerate(proc):
                for z in range(frames_per):
                    data_file[self.files.index(y[1])].append(f.read_frame(z + frames_per * ind))

        return data_file

    def swap_percent(self):

        if self.num_swaps == 1:
            raise ValueError("Swap was set to False")

        percent = []

        for r in self.info:
            swaps = 0
            for i in range(1, len(r)):
                if r[i][1] != r[i-1][1]:
                    swaps += 1
            ind1 = r[0][0].find('_R')
            ind2 = r[0][0].find('_', ind1 + 1)
            name = r[0][0][ind1 + 1:ind2]
            percent.append([name, swaps/(len(r)-1) * 100])
        return percent
