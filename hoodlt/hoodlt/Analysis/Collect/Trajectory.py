"""
:module: Trajectory
:platform: Unix, Windows
:synopsis: stores the reinitialization of a set of configurations from a gsd file 

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> April 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - units are now applied while reading in data from the gsd, not after
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Reworked class and moved it to the Analysis package
"""

import numpy as np
import gsd.hoomd
from hoodlt.Data.Modelconfigurations.Saver import load_config
from hoodlt.Analysis.Collect.ReInitHelper import reinit_frame


class Trajectory(object):
    """
    Stores the reinitialization of a set of configurations from a gsd file. Represents the entirety of a simulation.
    """

    def __init__(self, gsd_name, pickle_name, frames=None):
        """

        :param gsd_name: the name of the gsd file, without the .gsd extension
        :param pickle_name: the name of the pickle file, without the .pickle extension
        :param frames: list of frames to add to the trajectory, defaults to all frames in the gsd
        """

        sys = gsd.hoomd.open(gsd_name+'.gsd', mode='r')

        if frames is None:
            self.frames = list(range(len(sys)))
        else:
            self.frames = list(np.sort(frames))

        self.nframes = len(self.frames)
        self.timesteps = [sys[int(i)].configuration.step for i in self.frames]

        # build list of configurations
        conf_old = load_config(pickle_name)
        units = conf_old.ff_reader.get_units()
        self.configurations = []
        print("Reinitializing ...")
        for frame in self.frames:
            print("Frame:", frame)
            conf = reinit_frame(conf_old, sys, frame, units)
            self.configurations.append(conf)

    def dump(self, name):
        """
        Dumps the trajectory to a gsd file

        :param name: the name of the gsd file to dump to, without the .gsd extension
        :return: nothing
        """

        t = gsd.hoomd.open(name=name, mode='wb')

        t.extend([con.make_snapshot() for con in self.configurations])

    def append(self, trajectory2):
        """
        Appends all the data in the given trajectory to this trajectory

        :param trajectory2: trajectory to add, added to the end and frames and timestpes are re-calculated
        :return: void, this trajectory object will be changed
        """

        self.frames = self.frames + [f + self.frames[-1] for f in trajectory2.frames]
        self.nframes = self.nframes + trajectory2.nframes
        self.timesteps = self.timesteps + [t + self.timesteps[-1] for t in trajectory2.timesteps]
        self.configurations = self.configurations + trajectory2.configurations

    def remove_frames(self, frames):
        """
        Removes all the frames in the input list of frames from this trajectory object

        :param frames: list of the frames to be removed
        :return: void. Tracjectory Object is just changed
        """

        for frame in frames:
            if frame in self.frames:
                index = self.frames.index(frame)
                self.timesteps.remove(self.timesteps[index])
                self.frames.remove(frame)
                self.configurations.remove(self.configurations[index])


class ConfigurationTrajectory(Trajectory):
    """
    create a trajectory object out of configurations
    """
    def __init__(self, configs):
        """
        
        :param configs: list of FunctionalizedConfigurations for in the order they will be in in the trajectory
        """

        self.configurations = configs
        self.frames = list(range(len(self.configurations)))
        self.nframes = len(self.frames)
        self.timesteps = self.frames
