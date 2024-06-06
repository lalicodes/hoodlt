"""
:module: PlotSolventProperties
:platform: Unix, Windows
:synopsis: Defines the class used to plot Solvent Analysis Calculations and write data to files

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
"""

import numpy as np
import matplotlib.pyplot as plt
from hoodlt.Analysis.Analyze.AnalyzeSolvent import AnalyzeSolvent


class PlotSolventProperties:

    def __init__(self, comp):
        """

        :param comp: ComputeSolventProperties object
        """

        self.comp = comp

    def plot_local_density_hist(self, nbins=None, density_scale=1, liquid_dens=None, vapor_dens=None):
        """

        :param nbins: number of bins in the histogram
        :param density_scale: factor to scale all the density values by
        :param liquid_dens: optional experimental value of the liquid density
        :param vapor_dens: optional experimental value of the gas density
        :return: None
        """

        self._check_field("densities")

        plt.hist(density_scale * np.array(self.comp.densities).flatten(), bins=nbins)

        if liquid_dens is not None:
            plt.axvline(liquid_dens, color="green", alpha=.3, label="Experimental liquid")
        if vapor_dens is not None:
            plt.axvline(vapor_dens, color="red", alpha=.3, label="Experimental gas")

        plt.xlabel("Density")
        plt.ylabel("Number of Toluenes")
        plt.savefig("density_hist.png", bbox_inches="tight")
        plt.legend(loc="best")
        plt.show()

    def _check_field(self, field):
        """
        Makes sure the compute object has calculated some value

        :param field: the name of the field
        """

        if getattr(self.comp, field) is None:
            raise ValueError("Must compute the " + field + " in the ComputeSolventProperties object before this "
                                                           "operation can be done")

    def write_densities_to_file(self, file_name):
        """
        Writes all the densities calculated by the compute object to a text file

        :param file_name: name of the file, without the .txt extension, to write the data to
        :return: None
        """

        self._check_field("densities")

        dens = np.array(self.comp.densities)

        f = open(file_name+'.txt', mode='w')
        for i in range(len(dens[0])):
            for d in dens[:,i]:
                f.write("%s " % d)
            f.write("\n")

    def write_to_colored_xyz(self, frame, vapor_cut, liquid_cut, include_ncs=False):
        """
        Writes the positions of the solvents (and possibly ncs) to an xyz file where the solvents are colored by
        whether they are liquid, gas, or interface

        :param frame: frame to write to the xyz
        :param vapor_cut: vapor cutoff, in simulation units. Any solvent with local density less than the cutoff will
        not be considered part of the droplet. This number should be in simulation units
        :param liquid_cut: liquid cutoff, in simulation units. Any solvent with local density greater than the cutoff
        will be considered liquid. This number should be in simulation units
        :param include_ncs: Whether or not to include the nc positions in the dumped xyz file
        :return: None
        """

        densities, _ = self.comp.compute_local_density_by_frame(frame)
        config = self.comp.traj.configurations[self.comp.traj.frames.index(frame)]

        xyz_data = self._create_xyz_data(config, densities, vapor_cut, liquid_cut, include_ncs)

        # write the data to a file
        xyz = open("frame_"+str(frame)+"_liquid_vapor.xyz", 'w')
        xyz.write(str(len(xyz_data)))
        xyz.write("\n\n")
        for line in xyz_data:
            xyz.write("%s %s %s %s\n" %(line[0], line[1], line[2], line[3]))
        xyz.close()

    def _create_xyz_data(self, config, densities, vapor_cut, liquid_cut, include_ncs):
        """
        Creates the bulk of the xyz data to be written to a file

        :param config: the configuration to write to an xyz
        :param densities: a list of local densities of the solvents
        :param vapor_cut: vapor cutoff
        :param liquid_cut: liquid cutoff
        :param include_ncs: whether to inlcude the ncs in the xyz
        :return: a 2D list of xyz data
        """

        xyz_data = []

        # write the data for all the nanoparticles in the system
        if include_ncs:
            for nc in config.particles:
                self._add_nc_data(xyz_data, nc)

        # write the data for the solvents
        for i, dens in enumerate(densities):

            solv_pos = AnalyzeSolvent(config.solvent[i]).average_position()

            if dens > liquid_cut:
                self._add_entry_to_xyz(xyz_data, 'liquid', solv_pos)
            elif dens < vapor_cut:
                self._add_entry_to_xyz(xyz_data, 'vapor', solv_pos)
            else:
                self._add_entry_to_xyz(xyz_data, 'intermediate', solv_pos)

        return xyz_data

    def _add_nc_data(self, xyz_data, nc):
        """
        Adds the data for a given nc to the xyz_data array

        :param xyz_data: list which holds all the data to be written to the xyz
        :param nc: FuntionalizedParticle object
        :return: None
        """

        self._add_basic_entity_data_range(xyz_data, nc.core, nc.core.num_particles - nc.core.graft_num)
        for lig in nc.ligands:
            self._add_basic_entity_data_range(xyz_data, lig, lig.num_particles)

    def _add_basic_entity_data_range(self, xyz_data, basic_entity, size):
        """
        Adds the first size particles in the given basic entity to the xyz_data array

        :param xyz_data: an array holding data that will be written to the xyz
        :param basic_entity: BasicSystemEntity object
        :param size: the number of particles to add to the xyz
        :return: None
        """

        for i in range(size):
            self._add_entry_to_xyz(xyz_data, basic_entity.typeid[i], basic_entity.position[i])

    def _add_entry_to_xyz(self, xyz_data, type, position):
        """
        Adds a single entry to the xyz_data array

        :param xyz_data: the array that holds all the xyz data
        :param type: the type of the particle to add to the xyz
        :param position: the position of the particle to add to the xyz
        :return: None
        """

        xyz_data.append([type, position[0], position[1], position[2]])

    def plot_rdf_from_center_of_droplet(self, x_min, x_max, y_max, density_scale=1, distance_scale=1):
        """
        Makes a plot of the rdf data

        :param x_min: x min for the plot window
        :param x_max: x max for the plot window
        :param y_max: y max for the plot window
        :param density_scale: scale factor for the plotted densities
        :param distance_scale: scale factor for the plotted distances
        :return: None
        """

        self._check_field("rdf_data")

        (x, f_x) = self.comp.rdf_data

        x *= distance_scale
        f_x *= density_scale

        plt.plot(x, f_x)
        plt.xlim(x_min, x_max)
        plt.ylim(0, y_max)
        plt.xlabel(r"$R$")
        plt.ylabel(r"$\rho / \rho_{avg}$")
        plt.title(r"Density vs. distance from the center of the droplet")
        plt.legend(loc="upper right")
        plt.savefig("dens_from_center_of_drop_gibbs.png", bbox_inches="tight")
        plt.show()

    def plot_evaporation_rate(self, fit_func, time_scale=1):
        """
        Plots the fit of the number of gas particles along with the actual data for the number of gas particles

        :param fit_func: fitting function
        :param time_scale: factor by which to multiply all the times. If left to 1, the times will be range(len(num_gas))
        :return: None
        """

        self._check_field("num_gas")
        self._check_field("gas_fit_data")

        (popt, pcov) = self.comp.gas_fit_data

        times = range(len(self.comp.num_gas))
        times = [time_scale * t for t in times]

        fit_data = [fit_func(t, *popt) for t in times]

        # make plot with fit on top of data
        plt.errorbar(times, self.comp.num_gas, fmt='o', markersize=1, label="simulation data")
        plt.plot(times, fit_data, label="fit")
        plt.xlabel(r"$t$")
        plt.ylabel(r"$N_{gas}$")
        plt.legend(loc="lower right")
        plt.savefig("evap_rate_plot.png", bbox_inches="tight")
        plt.show()

    def write_gas_data_to_file(self, file_name):
        """
        Writes the data for the number of gas particles to a file

        :param file_name: name of the file to write to, without the .txt extension
        :return: None
        """

        self._check_field("num_gas")

        f = open(file_name+'.txt', 'w')
        for x in self.comp.num_gas:
            f.write("%s\n" % x)
        f.close()
