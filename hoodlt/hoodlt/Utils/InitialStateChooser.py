"""
:module: InitialStateChooser
:platform: Unix, Windows
:synopsis: finds and dumps states that are closest to rvalues needed for MPI runs 

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> June 2017
"""

from __future__ import division
import gsd.hoomd
import numpy as np
from numpy import linalg as la


class InitialStateChooser(object):
    """
    finds and dumps states that are closest to rvalues needed for MPI runs
    """

    def __init__(self, dump_file, rvals=None):
        """
        
        :param dump_file: string, the file from which the states are chosen
        :param rvals: list of numbers, if the file contains a serial run, it will store which rvals are contained
        :return: 
        """

        self.sys = gsd.hoomd.open(dump_file, mode='rb')
        self.filename = dump_file
        if rvals is None:
            self.rvals = None
        else:
            if len(self.sys) % len(rvals) != 0:
                raise ValueError("Illegal rvals. File contains " + str(len(self.sys)) + " but " + str(len(rvals))
                                 + "r values")
            rvals_per = len(self.sys) / len(rvals)
            self.rvals_per = rvals_per
            self.rvals = []
            for r in rvals:
                self.rvals += [r] * int(rvals_per)

    def choose_all(self, wanted_rs, conv_unit=1, cut=0):
        """
        
        :param wanted_rs:a list of values the same length and order as rvals
        :param conv_unit: unit to multiply distance units to be in wanted_rs units
        :param cut: a fraction between 0 and 1. Cuts off the first fraction of the frames to make sure an equiliibrated 
        state is found; defaults to 0
        :return: list of best found Rs in the same order as wanted rs and rvals. Also dumps the corresponding frames
        """

        if self.rvals is None:
            raise ValueError("Can not call this function without setting the rvals first")

        if len(wanted_rs) != len(list(set(self.rvals))):
            raise ValueError("wanted_rs must be as long as rvals: " + str(len(list(set(self.rvals)))))
        if cut >= 1 or cut < 0:
            raise ValueError("cut must be between greater than 0 and less than 1")

        found_rs = np.zeros(len(wanted_rs))

        for ind, r in enumerate(wanted_rs):
            close_frame = 0
            dist_apart = 100
            for ind2 in range(len(self.sys)):
                y = int(ind2/self.rvals_per)
                if y == ind and (ind2 % self.rvals_per) / self.rvals_per >= cut:
                    frame = self.sys.read_frame(ind2)
                    pos1 = frame.particles.position[0] + np.multiply(frame.particles.image[0],
                                                                     frame.configuration.box[:-3])
                    pos2 = frame.particles.position[1] + np.multiply(frame.particles.image[1],
                                                                     frame.configuration.box[:-3])
                    dist = la.norm(pos1 - pos2) * conv_unit
                    if abs(dist - r) < abs(dist_apart):
                        close_frame = ind2
                        dist_apart = dist - r
            found_rs[ind] = dist_apart + r
            ind1a = self.filename.find('_R')
            ind2a = self.filename.find('_', ind1a + 1)
            name1 = self.filename[:ind1a + 1]
            name2 = self.filename[ind2a:-16]
            t = type(self.rvals[close_frame])
            gsd.hoomd.create(name1+'R'+str(t(self.rvals[close_frame]))+name2+'_init.gsd',
                             snapshot=self.sys.read_frame(close_frame))

        return list(found_rs)

    def auto_choose(self, hist_data, conv_unit=1, cut=0):
        """
        
        :param hist_data: normal histogram data np.array((num rs, num points per r))
        :param conv_unit: conv_unit: unit to multiply distance units to be in wanted_rs units
        :param cut: a fraction between 0 and 1. Cuts off the first fraction of the frames to make sure an equiliibrated 
        state is found; defaults to 0
        :return: list of found Rs closest to the average of the window in the same order as wanted rs and rvals. 
        Also dumps the corresponding frames
        """

        if len(hist_data) != len(list(set(self.rvals))):
            raise ValueError("hist data must be as long as rvals: " + str(len(list(set(self.rvals)))))
        wanted_rs = [np.average(hist_data[i]) for i in range(len(hist_data))]
        print("searching for:")
        print(wanted_rs)

        return self.choose_all(wanted_rs, conv_unit=conv_unit, cut=cut)

    def choose_one(self, wanted_r, conv_unit=1):

        """
        searches the whole file for the state closest to a given r, returns the found r and dumps the file
        :param wanted_r: the r distance it wants to find
        :param conv_unit: unit to multiply distance units to be in wanted_rs units
        :return: the found r
        """

        close_frame = 0
        dist_apart = 100
        for ind in range(len(self.sys)):
            frame = self.sys.read_frame(ind)
            pos1 = frame.particles.position[0] + np.multiply(frame.particles.image[0], frame.configuration.box[:-3])
            pos2 = frame.particles.position[1] + np.multiply(frame.particles.image[1], frame.configuration.box[:-3])
            dist = la.norm(pos1 - pos2) * conv_unit
            if abs(dist - wanted_r) < abs(dist_apart):
                close_frame = ind
                dist_apart = dist - wanted_r

        ind1a = self.filename.find('_R')
        ind2a = self.filename.find('_', ind1a + 1)
        name1 = self.filename[:ind1a + 1]
        name2 = self.filename[ind2a:-16]
        gsd.hoomd.create(name1+'R'+str(wanted_r)+name2+'_init.gsd', snapshot=self.sys.read_frame(close_frame))
        return dist_apart + wanted_r
