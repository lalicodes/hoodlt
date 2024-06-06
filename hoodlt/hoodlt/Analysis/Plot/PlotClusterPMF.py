"""
:module: PlotClusterPMF
:Platform: Windows, Unix
:synopsis: Makes plots for a small cluster of nanoparticles

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - plots can be made in any energy unit
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - updated documentation
..                  - Moved old code to fit new analysis scheme
"""

import numpy as np
import matplotlib.pyplot as plt
from hoodlt.PhysChemValues.OPM import opm_val
from hoodlt.PhysChemValues.HydroCarbon import HydroCarbon
from hoodlt.PhysChemValues.Phys2lambxi import phys2lambxi


class PlotClusterPMF:
    """
    This class makes plots comparing the sum of pair potentials in small symmetric nc systems to the calculated pmf
    """
    def __init__(self, compute_pmf):
        """

        :param compute_pmf: object of type ComputeClusterPMF
        """

        self.avg_ctr_cage_dists = compute_pmf.compute_avg_ctr_cage_dists()
        self.free_energy = compute_pmf.compute_free_energy()
        self.sum_pair_potentials = compute_pmf.compute_sum_pair_potentials()
        (self.sppot_dists, self.symmetric_pair_potential) = compute_pmf.compute_symmetric_pair_potential()

        self.num_rvals = len(compute_pmf.rvals)
        self.num_cage_ncs = compute_pmf.cluster_obj.num_cage_ncs

    def plot_pmf(self, energy_scale=1):
        """
        Plots the computed pmf for the system

        :param energy_scale: factor by which to scale the values of energy in the plot, default is to plot in KbTs
        :return: None
        """

        free_energy = self.free_energy * energy_scale

        plt.plot(self.avg_ctr_cage_dists, free_energy)

        plt.savefig("pmf.png")
        plt.show()

    def write_pmf_plot_data(self, energy_scale=1):
        """
        Writes the data for the pmf to a file

        :param energy_scale: factor by which to scale the values of energy to write, default is to write in KbTs
        :return: None
        """

        fe = self.free_energy * energy_scale

        f = open("pmf_plot_data.txt", "w")
        for i in range(self.num_rvals):
            ind = self.num_rvals - 1 - i
            f.write("%s %s \n" % (self.avg_ctr_cage_dists[ind], fe[ind]))

    def plot_many_body_effects(self, x_min, x_max, num_core_atoms, core_diameter, graft_num, lig_chain_length, energy_scale=1):
        """
        Creates a plot comparing the computed pmf to the sum of pair potentials in the system

        :param x_min: minimum distance for the plot
        :param x_max: maximum distance for the plot
        :param num_core_atoms: number of atoms in the core of whichever nc was used in the simulation
        :param core_diameter: diameter of the core used in the simulation
        :param graft_num: number of chains grafted to the surface of the nc used in the simulation
        :param lig_chain_length: number of repeats in the ligand in the simulation. Ligand should be Hydrocarbon
        :param energy_scale: factor by which to scale the values of energy in the plot, default is to plot in KbTs
        :return: None
        """

        self._check_for_pair_data()

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=18)

        opm = self._compute_opm(lig_chain_length, core_diameter)

        # plot the data
        fe = self.free_energy * energy_scale
        spp = self.sum_pair_potentials * energy_scale

        plt.plot(self.avg_ctr_cage_dists, fe, "b-")
        plt.plot(self.avg_ctr_cage_dists, spp, "r-")
        plt.plot(self.avg_ctr_cage_dists, fe - spp, "g-")

        # add vertical lines at the minima of the curves
        plt.axvline(x=opm, color="black", linestyle="--")
        plt.axvline(x=self.avg_ctr_cage_dists[np.argmin(fe)], color="blue", linestyle="--")
        min_sum_pairs = np.argmin(spp)
        if min_sum_pairs != 0:
            plt.axvline(x=self.avg_ctr_cage_dists[min_sum_pairs], color="red", linestyle="--")
        plt.axhline(y=0, color="grey", alpha=.3)

        # add labels, legend, and title
        plt.tick_params(length=10, width=2.5, labelsize=16)
        plt.xlim(x_min, x_max)
        plt.xlabel(r"R (\AA)", size=18)
        plt.ylabel(r"$k_B T$", size=18)
        plt.legend(["$F_N(T,R)$", r"$\sum_{i \neq j}^N F_2(T, R_{ij})$", "$\Delta^{MB}(T,R)$", "OPM"],
                   loc="lower right")
        title_txt = "Au$_{" + str(num_core_atoms) + "}$(SC$_{" + str(lig_chain_length) + "}$)$_{" + str(
            graft_num) + "}$ \quad P$_" + str(self.num_cage_ncs) + "$"
        plt.title(title_txt, size=24)

        # save the figure
        plt.savefig("many_body_effects.png", bbox_inches="tight")
        plt.show()

    @staticmethod
    def _compute_opm(chain_length, diameter):
        """
        Computes the opm value of an nc with the given chain length and diameter

        :param chain_length: number of repeats on the ligand (it assumes Hydrocarbon)
        :param diameter: diameter of the core used in the simulation
        :return: the opm value
        """
        chain = HydroCarbon(chain_length)
        lamda, xi = phys2lambxi(10 * chain.max_length(), 1, .5 * diameter, 1)
        opm = opm_val(lamda, xi)
        opm = opm * diameter

        return opm

    def write_many_body_plot_data(self, energy_scale=1):
        """
        Writes the data for the many body plot to a file

        :param energy_scale: factor by which to scale the values of energy to write, default is to write in KbTs
        :return: None
        """

        self._check_for_pair_data()

        fe = self.free_energy * energy_scale
        spp = self.sum_pair_potentials * energy_scale

        # write plot data to files
        f = open("many_body_effects_plot_data.txt", "w")
        for i in range(self.num_rvals):
            ind = self.num_rvals - 1 - i
            f.write("%s %s %s %s\n" % (self.avg_ctr_cage_dists[ind], fe[ind], spp[ind], fe[ind] - spp[ind]))

    def _check_for_pair_data(self):
        """
        Checks if pair data is defined, if not, it throws and error

        :return: None
        """
        if self.sum_pair_potentials is None:
            raise ValueError("Must give a file with pair potential information")

    def plot_pair_pmf_comparison(self, energy_scale=1):
        """
        Makes a plot which compares the calculated sum of pair potentials to the perfectly symmetric sum of pair
        potentials

        :param energy_scale: factor by which to scale the values of energy in the plot, default is to plot in KbTs
        :return: None
        """

        self._check_for_pair_data()

        ppp = self.symmetric_pair_potential * energy_scale
        spp = self.sum_pair_potentials * energy_scale

        plt.plot(self.avg_ctr_cage_dists, spp, "r-")
        plt.plot(self.sppot_dists, ppp, "orange")

        plt.savefig("pair_potential_comparison.png")
        plt.show()

    def write_pair_pmf_comparison_plot_data(self, energy_scale=1):
        """
        Writes the data for the pair potential comparison plot to a file

        :param energy_scale: factor by which to scale the values of energy to write, default is to write in KbTs
        :return: None
        """
        self._check_for_pair_data()

        ppp = self.symmetric_pair_potential * energy_scale
        spp = self.sum_pair_potentials * energy_scale

        # write plot data to files
        f = open("pair_potential_comparison_plot_data.txt", "w")
        for i in range(len(self.sppot_dists)):
            if i < self.num_rvals:
                f.write("%s %s %s %s\n" % (self.sppot_dists[i], ppp[i], self.avg_ctr_cage_dists[self.num_rvals-1-i],
                                           spp[self.num_rvals-1-i]))
            else:
                f.write("%s %s\n" % (self.sppot_dists[i], ppp[i]))
