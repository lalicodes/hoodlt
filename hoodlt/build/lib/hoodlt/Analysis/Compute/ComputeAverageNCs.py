"""
:module: ComputeAverageNCs
:Platform: Windows, Unix
:synopsis: Computes averages of nanoparticles and configurations over the course of simulations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - class now explicitly takes a trajectory object
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Moved old code (written by Curt Waltmann) to fit new analysis scheme
..
..                Alex Travesset <trvsst@ameslab.gov>, August 2022
..                  - updated the code to make it usable given all changes made during 2022
"""

import copy as cp
import numpy as np

from hoodlt.Analysis.Compute.ComputeWithTrajectory import ComputeWithTrajectory
from hoodlt.Data.Modelconfigurations.FunctionalizedConfiguration import FunctionalizedConfiguration


class ComputeAverageNCs(ComputeWithTrajectory):
    """
    This class computes averages of NCs throughout the progression of simulations. It can make average NCs, average
    Configurations, or averages for each window in a simulation
    """
    def __init__(self, trajectory):
        """

        :param trajectory: a Trajectory object
        """

        super(ComputeAverageNCs, self).__init__(trajectory)

    def compute_average_particle(self, nano_index, frames=None, orientation_canonical=True):
        """
        Computes the average of a nanoparticle over the given number of frames

        :param nano_index: index of the particle to be averaged
        :param frames: list of the frames to be averaged, defaults to all that have been initialized
        :param orientation_canonical: True: positions in body frame. False: positions in average space frame
        :return: a FunctionalizedParticle object, centered at 0 with [1,0,0,0] orientation and average ligand positions
        """

        if nano_index >= len(self.traj.configurations[0].particles):
            raise IndexError("nano_index out of range")

        # check consistency of the frames
        frames = self._check_and_create_frames(frames)

        # get the box size
        box = self.traj.configurations[0].box
        # make a copy of the nc, and shift it to the origin
        part = cp.deepcopy(self.traj.configurations[0].particles[nano_index])
        part.shift(-1 * part.core.position[0])

        # define frame list
        go_through = [cp.deepcopy(self.traj.configurations[self.traj.frames.index(frame)]) for frame in frames]

        # define empty lists for averages
        ligand_positions = []
        linker_positions = []

        # obtain the lists
        for index, con in enumerate(go_through):
            con = cp.deepcopy(con)
            part = con.particles[nano_index]
            part.shift(-1 * part.core.position[0])
            if orientation_canonical:
                qrot = part.core.orientation[0]
                part.rotate_actual(qrot)
            linker_positions = [cp.deepcopy(lig.position[0]) for lig in part.ligands]
            for ind, lig in enumerate(part.ligands):
                lig.shift(-linker_positions[ind])
            ligand_positions.append([lig.position for lig in part.ligands])

        # compute average ligand position
        average_ligand = np.average(ligand_positions, axis=0)
        # restore positions
        for index, lig in enumerate(part.ligands):
            lig.position = average_ligand[index] + linker_positions[index]

        # return the nanoparticle as a functionalized configuration object centered at the origin
        vec = np.array([0.0, 0.0, 0.0])
        return FunctionalizedConfiguration([part], [vec], box)

    def compute_average_config(self, particle_indexes, frames=None, o_canon=False, centers=False):
        """
        Computes an average configuration of nanoparticles specified by the input list of indexes

        :param particle_indexes: indexes of the nanoparticles to include
        :param frames: list of frames to use
        :param o_canon: whether to return nanoparticles body frame
        :param centers: place the nanoparticles at the given centers
        :return: average configuration aver the frames
        """

        frames = self._check_and_create_frames(frames)
        box_l = self.traj.configurations[0].box
        mat = np.zeros([len(frames), len(particle_indexes), 3])
        center_defined = centers
        if not isinstance(centers, (bool,)):
            center_defined = True

        if not center_defined:
            go_through = [cp.deepcopy(self.traj.configurations[self.traj.frames.index(frame)]) for frame in frames]
            for index, con in enumerate(go_through):
                cn = cp.deepcopy(con)
                for ind, part in enumerate(particle_indexes):
                    mat[index, ind] = cn.particles[part].core.position[0]
            c_pos = np.average(mat, axis=0)
        else:
            c_pos = centers

        nanos = []
        for ind, part in enumerate(particle_indexes):
            nanos.append(self.compute_average_particle(part, frames=frames, orientation_canonical=o_canon).particles[0])

        return FunctionalizedConfiguration(nanos, c_pos, box_l)
