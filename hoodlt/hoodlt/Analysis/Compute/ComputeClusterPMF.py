"""
:module: ComputeClusterPMF
:Platform: Windows, Unix
:synopsis: Computes the Free Energy for a small cluster of nanoparticles

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Moved old code (written by Tommy Waltmann) to fit new analysis scheme
"""

import numpy as np
import scipy as sp
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Analysis.Collect.CollectHistData import CollectHistData
from hoodlt.Analysis.Collect.CollectGeneralData import CollectGeneralData
from hoodlt.Analysis.Collect.CollectForceFieldData import CollectForceFieldData


class ComputeClusterPMF:
    """
    This class computes the potential of mean force for small symmetric configurations of NCs and compares the pmf with
    the sum of pair free energies in the system.
    """

    def __init__(self, hist_files, rvals, cluster_obj, ff, bond_name_1, bond_name_2, pair_pmf_file=None,
                 relative_temp=1):
        """

        :param hist_files: a list of file names that store the histogram data from the simulation, without the .hist extension
        :param rvals: list of harmonic bond distances from the simulation. Input should be sorted descending
        :param cluster_obj: object that inherits from ClusterConfigDataAbs, and corresponds to the system that was run in the simulation
        :param units: units object. Units the simulation was run in
        :param ff: name of the forcefield used to run the simulation, without the '_forcefield.xlsx'
        :param bond_name_1: name of the center-cage bond used in the simulation
        :param bond_name_2: name of that cage-cage bond used in the simulation
        :param pair_pmf_file: file that contains the pair pmf
        :param relative_temp: temperature of the simulation divided by the default for the chosen unit system
        """

        ff_data = CollectForceFieldData(ff)
        hist_data = CollectHistData(hist_files, ff)
        other_data = CollectGeneralData(pair_pmf_file)

        if pair_pmf_file is not None:
            [self.pair_pmf_dists, self.pair_pmf_vals] = other_data.get_columns_from_text_file()

        units = ForceFieldReader(ff).get_units()  # TODO make sure calculations are done in simulation units
        self.avg_bond_lengths = hist_data.get_average_bond_dist_all_files()
        self.k1 = ff_data.get_bond_k(bond_name_1, units.energy_to_kbt(relative_temp))
        self.k2 = ff_data.get_bond_k(bond_name_2, units.energy_to_kbt(relative_temp))

        self.rvals = rvals
        self.cluster_obj = cluster_obj
        self.num_rvals = len(self.rvals)
        self.num_bonds = len(self.cluster_obj.bond_list)
        self.num_cage_ncs = cluster_obj.num_cage_ncs
        self.avg_ctr_cage_dists = self.compute_avg_ctr_cage_dists()
        self.avg_bond_displacement = self._compute_avg_bond_displacement()

    def compute_avg_ctr_cage_dists(self):
        """
        Computes the average center-cage distance for each simulation window

        :return: a 1D numpy array containing the average center-cage distance for each window of the simulation
        """

        plot_vals = np.zeros([self.num_rvals])
        for i in range(self.num_rvals):
            plot_vals[i] = np.average(self.avg_bond_lengths[i, :self.num_cage_ncs])

        return plot_vals

    def _compute_avg_bond_displacement(self):
        """
        Computes the average distance off the bond for each bond in each window of the simulation

        :return: A 2D numpy array with dimensions [number of bonds, number of windows]
        """

        delta_r = np.zeros([self.num_bonds, self.num_rvals])
        for rval in range(self.num_rvals):
            for bond in range(self.num_cage_ncs):
                delta_r[bond][rval] = self.avg_bond_lengths[rval][bond] - self.rvals[rval]

        return delta_r

    def compute_free_energy(self):
        """
        Computes the free energy, or pmf, of the system

        :return: a 1D numpy array whose values are the free energy of the system.
        """

        potential = np.zeros([self.num_rvals])

        for j in range(self.num_cage_ncs):
            for i in range(self.num_rvals):
                potential[i] += sp.trapz(self.avg_bond_displacement[j][i:self.num_rvals] *
                                         (self.k1 + (self.cluster_obj.bond_number_ratio*self.cluster_obj.bond_length_ratio**2) * self.k2),
                                         self.avg_ctr_cage_dists[i:self.num_rvals])

        # make the potential go to 0 as r gets very large
        potential -= potential[0]

        return potential

    def _pair_potential(self, avg_dist):
        """
        Computes the pair potential for a given list of distances

        :param avg_dist: a list of average distances for a bond in the system
        :return: a 1D numpy array which contains the values of the pair potential for a given bond in the system as the simulation progresses
        """

        return np.interp(avg_dist, self.pair_pmf_dists, self.pair_pmf_vals)

    def compute_sum_pair_potentials(self):
        """
        Computes the sum of all the pair free energies in the system

        :return: a 1D numpy array containing the values of the total pairwise free energy in the system
        """

        self._check_for_pair_data()

        potential = np.zeros([self.num_rvals])
        print(self.num_bonds)
        for bond in range(self.num_bonds):
            potential += self._pair_potential(self.avg_bond_lengths[:,bond])

        return potential

    def compute_symmetric_pair_potential(self):
        """
        Computes the what the pair potentials in the system would be if the system had stayed perfectly symmetric
        throughout the course of the simulation. Note: this requires no knowledge of the simulation to compute.

        :return: a tuple (domain , potential) where both are 1D arrays and domain is the x-values on which the values of potential should be plotted
        """

        sppot = []  # symmetric pair potentials
        domain = list(range(10 * self.rvals[-1], 10 * self.rvals[0] + 1))
        for i in range(len(domain)):
            domain[i] = domain[i] / 10

        for i in range(len(domain)):
            sppot.append(self.num_cage_ncs * self._pair_potential(domain[i]) +
                         (self.num_bonds - self.num_cage_ncs) * self._pair_potential(self.cluster_obj.bond_length_ratio * domain[i]))

        return np.array(domain), np.array(sppot)

    def _check_for_pair_data(self):
        """
        Checks if the pair potentials are defined or not. If not, it throws an error. This is a helper which is only
        called before a calculation involving pair potentials is done.

        :return: None
        """

        if self.pair_pmf_dists is None:
            raise ValueError("Must give a file with pair potential information")
