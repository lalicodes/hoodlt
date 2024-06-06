"""
:module: PlotWhamPotential
:platform: Unix, Windows
:synopsis: Defines the class used to plot PMF results

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - all energy plots can now be plotted in any energy unit
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - updated documentation
..                  - Moved old code (written by Alex Travesset) to fit new analysis scheme
..                Alex Travesset <trvsst@ameslab.gov>, February 2022
..                  - Added function to plot approximate pmf
..                  - updated documentation
..                Jonas Hallstrom <jonasleo@hotmail.com>, June 2022
..                  - Updated plot formatting in the plot_histograms function
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as line


class PlotWhamPotential:
    """
    This class plots wham potentials and other thermodynamic functions
    """
    def __init__(self, wham_class):
        """

        :param wham_class: object of type ComputeWhamPotential
        """

        self.wham = wham_class
        self.wham2 = self.wham.energies
        if self.wham.solution is None:
            print(" WARNING. The wham solution must be computed to be graphed. Can not graph wham PMF, "
                  "only approximate PMF.")

        self.data = self.wham.hist_data

    def plot_pmf(self, energy_scale=1):
        """
        Plots the pmf as computed by the WHAM method

        :param energy_scale: factor by which to scale the values of energy in the plot, default is to plot in KbTs
        :return: None
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('R')

        if energy_scale == 1:
            ax.set_ylabel('PMF (k$_{B}$T)')
        else:
            ax.set_ylabel('PMF')

        ax.set_title('PMF')
        pf = self.wham.pmf() * energy_scale
        val_min = np.amin(pf)
        ax.plot(self.wham.r, (pf-val_min))
        plt.savefig('PMF.png')
        plt.show()

    def plot_pmf_approx(self, energy_scale=1):
        """
        Plots the approximate pmf providing a cross-check for WHAM or useful when it cannot be solved

        :param energy_scale: factor by which to scale the values of energy in the plot, default is to plot in KbTs
        :return:
        """

        x, y = self.wham.PMF_approx()
        y *= energy_scale
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('R')

        if energy_scale == 1:
            ax.set_ylabel('PMF (k$_{B}$T)')
        else:
            ax.set_ylabel('PMF')

        ax.set_title('PMF')
        val_min = np.amin(y)
        ax.plot(x, (y-val_min))
        plt.savefig('PMF_approx.png')
        plt.show()

    def plot_prob_window(self, window_index='all', separate_windows=True):
        """
        Compute the probability that distance :math:'R' is measured in window i

        :param window_index: states which windows to plot if a numpy array
        :param separate_windows: Plot on separate windows
        """

        prob = self.wham.prob()
        num_windows, num_points = prob.shape

        if isinstance(window_index, np.ndarray):
            num_list = window_index
        else:
            num_list = range(num_windows)

        if separate_windows:
            for num in num_list:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_xlabel('R')
                win = 'window no' + str(num)
                plt.plot(self.wham.r, prob[num], label=win)
                plt.legend()
                f_name = 'prob ' + win + '.png'
                plt.savefig(f_name)
                plt.show()

        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('R')
            for num in num_list:
                win = 'window no' + str(num)
                plt.plot(self.wham.r, prob[num], label=win)
            plt.legend()
            f_name = 'prob' + '.png'
            plt.savefig(f_name)
            plt.show()

    def compare_bias_distributions(self):
        """
        Makes many plots comparing the actual distributions of histogram data to the computed distributions

        :return: graphs of computed vs. measured biased distributions for each window
        """

        for ind in range(len(self.wham.solution)):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('R')
            plt.plot(self.wham.r, self.wham.hist_simulation(ind), label='measured')
            plt.plot(self.wham.r, self.wham.hist_predicted(ind), label='computed')
            ax.set_title('computed/measured bin %d' % (ind + 1))
            plt.legend()
            plt.show()

    def plot_thermo_quantities(self):
        """
        Plots the free energy, entropy, and internal energy of the system

        :return:
        """

        free_energies = self.wham.free_energies()
        internal_energy = self.wham.internal_energy()
        entropies = self.wham.entropies()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('R')
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=12)

        ax.set_ylabel('Energy (k$_{B}$T)')

        ax.set_ylabel('Energy (k$_{B}$T)')
        ax.plot(self.wham.r, free_energies, label='Free Energy(PMF), F')
        ax.plot(self.wham.r, internal_energy - internal_energy[-1],
                label='Internal Energy, U', color='orange')
        ax.plot(self.wham.r, entropies - entropies[-1], label='Entropy, TS', color='green')
        plt.axhline(y=0, color="grey", alpha=.3)
        plt.legend(loc=4)
        plt.savefig('thermo.png')
        plt.show()

    def plot_2d_histogram(self, all=False):
        """
        plots the 2d histogram of r and internal energy

        :param all: when true every individual window is shown, when false just the combined window
        :return:
        """

        if self.wham2 is None:
            raise ValueError("Must call add2d before using this function, plot_2d_histogram")

        if len(self.wham.hist_2d.shape) != 3:
            raise ValueError('Must be initialized with 2d wham object to call this function')
        if all:
            for i in range(self.wham.num_windows):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title('2d histogram: Window ' + str(i))
                X, Y = np.meshgrid(self.wham.u, self.wham.r)
                ax.pcolormesh(X, Y, self.wham.hist_2d[i])
                plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('2d histogram: All Windows')
        X, Y = np.meshgrid(self.wham.u, self.wham.r)
        ax.pcolormesh(X, Y, np.sum(self.wham.hist_2d, axis=0))
        plt.savefig('2d.png')
        plt.show()

    def write_internal_energy(self, filename='internal_energy.txt', energy_scale=1):
        """
        writes the data for the calculated internal energy as a function of distance to be used in final figure

        :param filename: name of the file to write the data to
        :param energy_scale: factor by which to scale the values of energy, default is to write the data in KbT's
        """

        p = self.wham.internal_energy() * energy_scale
        fid = open(filename, 'w')
        for x in range(len(self.wham.r)):
            fid.write(str(self.wham.r[x]) + " " + str(p[x]))
            fid.write('\n')

    def write_entropy(self, filename='entropy.txt', energy_scale=1):
        """
        writes the data for the calculated entropy as a function of distance to be used in final figure generation

        :param filename: name of the file to write the data to
        :param energy_scale: factor to scale the values of energy, default is to write the data in KbT's
        """

        p = self.wham.entropies() * energy_scale
        fid = open(filename, 'w')
        for x in range(len(self.wham.r)):
            fid.write(str(self.wham.r[x]) + " " + str(p[x]))
            fid.write('\n')

    def write_pmf(self, filename='PMF.txt', energy_scale=1):
        """
        writes the data for the calculated PMF as a function of distance to be used in final figure generation

        :param filename: name of the file to write the data to
        :param energy_scale: factor to scale the values of energy, default is to write the data in KbT's
        """

        p = self.wham.free_energies() * energy_scale
        fid = open(filename, 'w')
        for x in range(len(self.wham.r)):
            fid.write(str(self.wham.r[x]) + " " + str(p[x]))
            fid.write('\n')

    def write_pmf_approx(self, filename='PMF_approx.txt', energy_scale=1):
        """
        writes the data for the calculated PMF using the approximate pmf calculation

        :param filename: name of the file to write data to
        :param energy_scale: factor to scale the values of energy, default is to write the data in kBT's
        """

        x, y = self.wham.PMF_approx()*energy_scale
        const = y[-1]
        fid = open(filename, 'w')
        for ind in range(self.wham.num_windows):
            z = self.wham.rvals[ind]
            # it plots r_0=original value for spring term, r, r_0-r, pmf
            fid.write('%1.7lf %1.7lf %1.7lf %1.7lf \n' % (z, x[ind], z-x[ind], y[ind]-const))

    def plot_histograms(self, r_values, binwidth=0.2):
        """
        creates the plot of histograms as a function of distance for all the bond lengths
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=12)

        colors = ['g', 'brown', 'b', 'r', 'cyan', 'magenta'] * len(self.data)

        maxd = 0
        mind = 1000
        for i in range(len(self.data)):
            newdat = self.data[i]
            maxd = max([max(newdat), maxd])
            mind = min([min(newdat), mind])
            plt.hist(newdat, bins=np.arange(min(newdat), max(newdat) + binwidth, binwidth),
                     color=colors[i], alpha=0.5, label="R=%d" % r_values[i])
        for i in range(len(self.data)):
            ax.scatter([r_values[i]], [-10], c='k', marker='o', s=300)
            ax.scatter([r_values[i]], [-10], c=colors[i], marker='o', s=250)

        ax.plot(r_values, [0] * len(r_values), 'k')
        l = line.Line2D([mind - 2, maxd + 2], [0, 0], color='k', linewidth=2)
        ax.add_line(l)
        ax.set_xlim([mind - 2, maxd + 2])
        ax.set_ylim([-20, plt.ylim()[1]])
        ax.set_xlabel(r'$\rho$')
        ax.set_ylabel('counts')
        ax.set_title('PMF Histograms')

        leg = ax.legend(bbox_to_anchor=(0., 1.06, 1., 0.204), loc='lower left', ncol=8, mode="expand", fontsize="xx-small")
        plt.savefig("histograms.png", bbox_inches='tight', dpi=200)
        plt.show()
