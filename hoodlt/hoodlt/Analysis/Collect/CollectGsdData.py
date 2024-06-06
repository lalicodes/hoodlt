"""
:module: CollectGsdData
:Platform: Windows, Unix
:synopsis: Collects data from a gsd file and applies units

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - eliminated explicit references to units classes
"""

from hoodlt.Analysis.Collect.Trajectory import Trajectory


class CollectGsdData:
    """
    This class collects data from gsd files and can make trajectory objects from gsd files. All objects/data returned
    will be in simulation units (not dimensionless units)
    """
    def __init__(self, gsd_name):
        """

        :param gsd_name: name of the gsd file, without the '.gsd' extension
        """

        self.gsd_name = gsd_name

    def create_trajectory(self, pickle_name, frames):
        """
        Makes a trajectory object which has all the chosen frames

        :param pickle_name: name of the pickle file created before the simulation was run, at the same time the initial configuration was first dumped. Without the .pickle extension
        :param frames: a list of frames to make into a trajectory object
        :return: a Trajectory object
        """

        return Trajectory(self.gsd_name, pickle_name, frames=frames)
