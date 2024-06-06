"""
:module: ComputeSolventProperties
:Platform: Windows, Unix
:synopsis: Computes Properties of the Solvents in a simulation

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added many different calculations to the class
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Moved old code to fit new analysis scheme
"""

import numpy as np
import numpy.linalg as la
from scipy.optimize import curve_fit
from hoodlt.Analysis.Collect.CollectGeneralData import CollectGeneralData
from hoodlt.Analysis.Compute.ComputeWithTrajectory import ComputeWithTrajectory
from hoodlt.Analysis.Analyze.AnalyzeSolvent import AnalyzeSolvent


class ComputeSolventProperties(ComputeWithTrajectory):
    """
    this class uses trajectory objects to compute properties of solvent
    """
    def __init__(self, trajectory):
        """

        :param trajectory: a Trajectory object
        """

        super(ComputeSolventProperties, self).__init__(trajectory)

        self._adjust_positions()

        self.densities = None  # this gets initialized once the local densities are computed
        self.volumes = None # this gets initialized once the local densities are computed
        self.num_gas = None  # this gets initiailized once compute_num_gas is called
        self.rdf_data = None  # this gets initialized once the rdf has been computed
        self.gas_fit_data = None

    def compute_local_density(self, frames=None):
        """
        Uses the voronoi construction to compute the local density of the solvent in this configuration. The volumes
        used to compute the ouput histogram will be an accumulation of the volumes in every frame of this trajectory
        object. Make sure you have freud installed with version 0.8.2 or later

        :param frames: list of frames to do the local density calculation over
        :return: a 2D numpy array which has all the densities for each frame as calculated by the voronoi construction
        """

        list_frames = self._check_and_create_frames(frames)

        print("Computing Local Densities ...")
        densities_all_frames = []
        volumes_all_frames = []
        for frame in list_frames:
            d, v = self.compute_local_density_by_frame(frame)
            densities_all_frames.append(d)
            volumes_all_frames.append(v)

        self.densities = np.array(densities_all_frames)
        self.volumes = np.array(volumes_all_frames)

        return self.densities

    def compute_local_density_by_frame(self, frame):
        """
        Computes the local density via voronoi construction for a single frame of the trajectory

        :param frame: the frame number
        :return: a list of local densities of the solvent molecules in the configuration, and the associated volumes
        """

        print("Frame: ", frame)
        config = self.traj.configurations[self.traj.frames.index(frame)]
        center_pos = []
        for solv in config.solvent:
            center_pos.append(AnalyzeSolvent(solv).average_position())
        b = box.Box(Lx=config.box.Lx, Ly=config.box.Ly, Lz=config.box.Lz)
        max_dim = np.max([config.box.Lx, config.box.Ly, config.box.Lz])
        vor = voronoi.Voronoi(b, max_dim / 2)
        vor.compute(np.array(center_pos), b, max_dim / 2)
        vor.computeVolumes()
        volumes = vor.getVolumes()
        densities = [config.solvent[j].get_mass() / vol for j, vol in enumerate(volumes)]

        return densities, volumes

    def _adjust_positions(self):
        """
        helper function for local_solvent_density(). This CHANGES POSITIONS of solvent molecules which are across the
        box

        :return: None
        """

        for config in self.traj.configurations:
            for solv in config.solvent:
                if AnalyzeSolvent(solv).max_radius() > config.box.Lx / 2:  # determine if the molecule positions are split across the box
                    ctr = AnalyzeSolvent(solv).average_position()
                    for pos in solv.position:
                        if pos[0] - ctr[0] > config.box.Lx/2:
                            pos[0] -= config.box.Lx
                        if pos[0] - ctr[0] < -config.box.Lx/2:
                            pos[0] += config.box.Lx
                        if pos[1] - ctr[1] > config.box.Ly/2:
                            pos[1] -= config.box.Ly
                        if pos[1] - ctr[1] < -config.box.Ly/2:
                            pos[1] += config.box.Ly
                        if pos[2] - ctr[2] > config.box.Lz/2:
                            pos[2] -= config.box.Lz
                        if pos[2] - ctr[2] < -config.box.Lz/2:
                            pos[2] += config.box.Lz

    def compute_rdf_from_center_of_droplet(self, liquid_cut, frames=None, x_max=60, delta_r=.01):
        """
        Computes the rdf as a function of distance from the center of the solvent droplet. You must call
        compute_local_density before this calculation is done. The list of frames provided as an argument should match
        the frames (in the correct order) used to compute the local densities.

        :param liquid_cut: liquid cutoff, in simulation units. Any solvent with local density less than the cutoff will not be considered for computing the center of the droplet. This number should be in simulation units
        :param file_name: name of the file with local densities, if you have already computed the local density and don't want to do it again
        :param x_max: maximum distance from the droplet center to compute the rdf for
        :param delta_r: width of the bins for the rdf histogram calculation
        :return: a tuple (x, y) which is the x and y values for the rdf calculation
        """

        list_frames = self._check_and_create_frames(frames)

        if self.densities is None:
            raise ValueError("Must call compute_local_density() before this calculation can be done")

        rdf = density.RDF(rmax=x_max, dr=delta_r)
        for ind, frame in enumerate(list_frames):

            # get the frame
            frame_data = self.traj.configurations[self.traj.frames.index(frame)]

            # build a list of positions to represent the solvent positions
            positions = [AnalyzeSolvent(solv).average_position() for solv in frame_data.solvent]

            # compute the center of the droplet
            pos_center_drop = self._compute_center_droplet_position(self.densities[ind], np.array(positions), liquid_cut)

            # find closest solvent to the center of the droplet
            pos_closest = self._find_closest_to_droplet_center(positions, pos_center_drop)

            # accumulate rdf using the closest position to the droplet center as the reference point
            b = self.traj.configurations[ind].box
            sbox = box.Box(Lx=b.Lx, Ly=b.Ly, Lz=b.Lz)
            rdf.accumulate(sbox, np.array([pos_closest]), np.array(positions))

            self.rdf_data = (rdf.getR(), rdf.getRDF())

        return self.rdf_data

    def _find_closest_to_droplet_center(self, positions, pos_center_drop):
        """
        Finds the position in the list of positions that is closest to pos_center_drop

        :param positions: list of positions
        :param pos_center_drop: center of the droplet
        :return: the closest position in the list of positions to the center of the droplet
        """

        dist_min = 10000000.0  # some arbitrary large number
        pos_closest = None
        for pos in positions:
            dist_this = la.norm(pos - pos_center_drop)
            if dist_this < dist_min:
                dist_min = dist_this
                pos_closest = pos

        return pos_closest

    def _compute_center_droplet_position(self, densities, positions, liquid_cut):
        """
        Helper to compute the positions of the center of the droplet for a given frame

        :param densities: list of densities for the given frame
        :param positions: list of positions for the given frame
        :param liquid_cut: liquid cutoff, any solvent which has local density higher than the cutoff is considered a part of the droplet. This number should be in simulation units
        :return:
        """

        # compute the position of the center of the droplet
        pos_center_drop = np.zeros(3)
        num_liq = 0
        for i, dens in enumerate(densities):
            if dens >= liquid_cut:
                num_liq += 1
                pos_center_drop = np.add(pos_center_drop, positions[i])

        return np.divide(pos_center_drop, num_liq)

    def compute_droplet_radius(self, vapor_cut, liquid_cut, frames=None):
        """
        Computes the droplet radius considering only liquid particles, and the droplet radius considering both liquid
        and intermediate particles. Scaling to be consistent with the Gibbs convention is left to the user.

        :param vapor_cut: vapor cutoff, in simulation units. Any solvent with local density less than the cutoff will not be considered part of the droplet. This number should be in simulation units
        :param liquid_cut: liquid cutoff, in simulation units. Any solvent with local density greater than the cutoff will be considered liquid. This number should be in simulation units
        :param frames: list of frames over which to calculate the droplet radius, defaults to all possible that were reinitialized
        :return: the average radius over the list of frames considering only liquid particles, and the average radius considering both liquid and interface particles
        """

        list_frames = self._check_and_create_frames(frames)

        drop_rads_liq = []
        drop_rads_int = []

        for ind, frame_nbr in enumerate(list_frames):

            if self.densities is None:
                densities, volumes = self.compute_local_density_by_frame(frame_nbr)
            else:
                densities = self.densities[list_frames.index(frame_nbr)]
                volumes = self.volumes[list_frames.index(frame_nbr)]

            drop_vol_liq = 0
            drop_vol_int = 0  # droplet volume including intermediate density solvents

            # get volume of droplets
            for i, dens in enumerate(densities):
                if (dens >= liquid_cut):
                    drop_vol_liq += volumes[i]
                if (dens >= vapor_cut):
                    drop_vol_int += volumes[i]

            drop_rads_liq.append((3 / 4 / np.pi * drop_vol_liq) ** (1 / 3))
            drop_rads_int.append((3 / 4 / np.pi * drop_vol_int) ** (1 / 3))

        return np.average(drop_rads_liq), np.average(drop_rads_int)

    def compute_solvent_fraction_in_cone(self, ocm_angle, frames=None):
        """
        Computes the the fraction of solvents which lie within a cone defined by the given angle. This is intended to
        be used on configurations of 2 nanoparticles restricted to move on the x-axis

        :param ocm_angle: angle to use for the definition of the cone
        :param frames: list of frames over which to calculate an average, defaults to all possible that were reinitialized
        :return:
        """

        list_frames = self._check_and_create_frames(frames)

        list_num_cone_solvs = []
        for frame in list_frames:

            # get the configuration
            conf = self.traj.configurations[self.traj.frames.index(frame)]

            # get the positions of the core centers
            pos_core = conf.particles[0].core.position[0]
            pos_core_2 = conf.particles[1].core.position[0]

            # determine how many are within the cone
            num_cone_solvs = 0
            for solv in conf.solvent:
                if self._is_within_cone(AnalyzeSolvent(solv).average_position(), pos_core, pos_core_2, ocm_angle):
                    num_cone_solvs += 1

            list_num_cone_solvs.append(num_cone_solvs)

        return np.average(list_num_cone_solvs) / len(self.traj.configurations[0].solvent)

    @staticmethod
    def _is_within_cone(pos_particle, pos_center_0, pos_center_1, theta):
        """
        Determines whether or not the given particle is within the cone defined by the 3 other parameters

        :param pos_particle: position of the particle
        :param pos_center_0: position of one of the core centers
        :param pos_center_1: position of the other core center
        :param theta: the angle which defines the cone
        :return: Whether the particle is within the cone defined by the angle and the two positions
        """

        if (pos_particle[0] < pos_center_0[0] and pos_particle[0] < pos_center_1[0]) or (
                pos_particle[0] > pos_center_0[0] and pos_particle[0] > pos_center_1[0]):
            return False

        if la.norm(pos_particle - pos_center_0) < la.norm(pos_particle - pos_center_1):
            closer_core_pos = pos_center_0
        else:
            closer_core_pos = pos_center_1

        # assume ncs are constrained to x axis
        calc_theta = np.arctan(np.sqrt(pos_particle[1] ** 2 + pos_particle[2] ** 2) /
                               np.abs(pos_particle[0] - closer_core_pos[0]))

        return calc_theta <= theta

    def compute_num_gas(self, vapor_cut):
        """
        Compute the number of gas particles in the system over a range of frames

        :param vapor_cut: cutoff for determining which particles are gas and which are not. Any particle with local density lower than the vapor cut will be considered gas particles. This number should be in simulation units
        :return: a list of the number of gas particles in each frame considered for the calculation
        """

        if self.densities is None:
            raise ValueError("Must call compute_local_density() before this calculation can be done")

        list_num_gas = []
        for densities in self.densities:

            num_gas = 0
            for d in densities:
                if d < vapor_cut:
                    num_gas += 1
            list_num_gas.append(num_gas)

        self.num_gas = list_num_gas

        return list_num_gas

    def compute_relaxation_parameters(self, fit_func, x_scale=1, frames=None, file_name=None):
        """
        Computes a curve fit to the data from the compute_num_gas() function. You must call compute_num_gas() before
        calling this method, or supply data from a single column text file.

        :param fit_func: function to use for fitting the gas data
        :param x_scale: scaling factor for the x-axis values, if left to default, the x-axis values will be range(len(frames))
        :param frames: list of frames over which the gas data was calculated, defualts to all possible
        :param file_name: name of the file used for supplying gas data, should be a single column
        :return: the parameters for the fitting function as returned by scipy.optimize.curve_fit
        """

        list_frames = self._check_and_create_frames(frames)

        if file_name is not None:
            self.num_gas = CollectGeneralData(file_name).get_columns_from_text_file()[0]
        elif self.num_gas is None:
            raise ValueError("Must call compute_num_gas() or supply data from a text file before this calculation can "
                             "be done")

        x_vals = [x_scale*n for n in range(len(list_frames))]

        self.gas_fit_data = curve_fit(fit_func, x_vals, self.num_gas)

        return self.gas_fit_data
