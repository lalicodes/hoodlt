"""
:module: ComputeWithTrajectory
:Platform: Windows, Unix
:synopsis: Computes simulation quantities by analyzing data in gsd files

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Class now explicity takes a trajectory object as an argument
"""


class ComputeWithTrajectory:
    """
    Abstract class that will serve as a super class for every class that does calculations on trajectory objects.
    """
    def __init__(self, trajectory):
        """

        :param trajectory: a Trajectory object
        """

        self.traj = trajectory

    def _check_and_create_frames(self, frames):
        """
        checks if the list of frames given by the user is valid

        :param frames: a list of frames
        :return: a list of frames to be used for various types of calculations
        """

        if frames is None:
            frames = self.traj.frames

        if not set(frames).issubset(set(self.traj.frames)):
            print(frames)
            print(self.traj.frames)
            raise ValueError("All frames to be average must have been initialized")

        return frames
