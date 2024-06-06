"""
:module: PlotLigandProperties
:platform: Unix, Windows
:synopsis: Defines the class used to plot Ligand calculations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - updated documentation
..                  - Moved old code (written by Curt Waltmann) to fit new analysis scheme
"""

import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from hoodlt.Analysis.Analyze.AnalyzeCore import AnalyzeCore
from hoodlt.Analysis.Analyze.AnalyzeLigand import AnalyzeLigand
from hoodlt.Analysis.Analyze.AnalyzeNanoparticle import AnalyzeNanoparticle as AnalyzeNC


class PlotLigandProperties:
    """
    This class makes plots of calculations done on the ligands
    """
    def __init__(self, compute, rvals=None):
        """

        :param compute: object of type ComputeLigandProperties
        :param rvals: optional list of rvals needed for some plots
        """
        self.compute = compute
        self.rvals = rvals

        if rvals is not None:
            if self.compute.traj.nframes % len(rvals) != 0:
                raise ValueError(str(len(rvals)) + " r values is not allowed for " + str(self.compute.traj.nframes)
                                 + "frames.")
            frames_per = self.compute.traj.nframes / len(rvals)
            self.rvals = self.rvals * int(frames_per)
            self.rvals.sort(key=float, reverse=True)
            self.rvals = [r for ind, r in enumerate(self.rvals) if ind in self.compute.traj.frames]

    def trans_vs_timestep(self, particle_indexes, ns=[None]):
        """
        Plots the percentage of trans dihedrals vs the timestep of the simulation

        :param ns: list of dihedral indexes to get graph the percent trans for, defaults to calculating them all as one
        :param particle_indexes: list of indexes of the nano particle to do the calculation for
        :return: plots the percent trans on the particle vs the timestep from the trajectory
        """

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        var = 'Percent Trans'

        ax1.set_title(var)
        ax1.set_xlabel('Timestep')
        ax1.set_ylabel(var)

        for i in particle_indexes:
            for n in ns:
                y = [AnalyzeNC(conf.particles[i]).percent_trans(n=n) for conf in self.compute.traj.configurations]
                label = str(i) + 'th particle'
                if n is not None:
                    label += ", " + str(n) + "th dihedral"
                ax1.plot(self.compute.traj.timesteps, y, label=label)
        plt.legend()
        plt.savefig('trans_v_timestep.png')
        plt.show()

    def distance_v_timestep(self, particle_indexes, ns):
        """
        Plots the average distance between the center of mass of the core and the particle of index n on the ligand vs
        the timestep of the simulation

        :param particle_indexes: list of indexes of the nano particle to do the calculation for
        :param ns: list of chain posiitions to do the calculation for
        :return: plots the average distance the nth particle fro the ligands on the particles vs the timestep of the trajectory dump step
        """
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        var = 'Distance to particle'

        ax1.set_title(var)
        ax1.set_xlabel('Timestep')
        ax1.set_ylabel(var)

        for i in particle_indexes:
            for n in ns:
                y = [AnalyzeNC(conf.particles[i]).avg_distance_nth_particle(n)
                     for conf in self.compute.traj.configurations]
                ax1.plot(self.compute.traj.timesteps, y, label=str(i) + 'th nanoparticle ' + str(n) + ' in chain')

        plt.legend()
        plt.savefig('distance_v_timestep.png')
        plt.show()

    def inter_angle(self, particle_indexes, frames):
        """
        the intermolecular angle is the angle between the chain vectors for every combination of chains on the particle.
        The chain vector is simply the vector between the first and last particles in the chain

        :param particle_indexes: list of indexes of the nano particle to do the calculation for
        :param frames: list of frames to do the calculation for
        :return: plots the intermolecular angle distributions for given frames and particles
        """

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        var = 'Intermolecular Angle Distribution'
        ax1.set_title(var)
        ax1.set_xlabel('Theta')
        ax1.set_ylabel('Intensity')

        colors = ['g', 'b', 'r', 'cyan'] * len(frames) * len(particle_indexes)

        for pi, i in enumerate(particle_indexes):
            for px, x in enumerate(frames):
                data = AnalyzeNC(self.compute.traj.configurations[self.compute.traj.frames.index(x)]).intermolecular_angles()
                hist, bin_edges = np.histogram(data, bins=np.arange(0, 190, 10))
                plt.plot(np.arange(5, 185, 10), hist / sum(hist), color=colors[px + pi], label=str(i) + "th particle, "
                                                                                            "frame " + str(x))
        plt.legend()
        plt.show()

    def trans_v_r(self, particle_indexes, ns=[None]):
        """
        trans is defined as a dihedral angle between + and - pi/3

        :param ns: list of dihedral indexes to get graph the percent trans for, defaults to calculating and graphing them all as one
        :param particle_indexes: list of indexes of the nano particle to do the calculation for
        :return: plots the average distance the nth particle fro the ligands on the particles for different windows
        """
        if self.rvals is None:
            raise AttributeError("R_values were never established in constructor")

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        var = 'Percent Trans'

        ax1.set_title(var)
        ax1.set_xlabel('R')
        ax1.set_ylabel(var)

        set_r = list(set(self.rvals))
        set_r.sort(key=float, reverse=True)
        for i in particle_indexes:
            for n in ns:
                trans = []
                for r in set_r:
                    y = [AnalyzeNC(self.compute.traj.configurations[t].particles[i]).percent_trans(n=n)
                         for t in range(len(self.compute.traj.configurations)) if self.rvals[t] == r]
                    trans.append(np.average(y))
                label = str(i) + 'th particle'
                if n is not None:
                    label += ", " + str(n) + "th dihedral"
                ax1.plot(set_r, trans, label=label)
        plt.legend()
        plt.show()

    def distance_v_r(self, particle_indexes, ns):
        """

        :param particle_indexes: list of indexes of the nano particle to do the calculation for
        :param ns: list of chain posiitions to do the calculation for
        :return: plots the average distance the nth particle fro the ligands on the particles
        """
        if self.rvals is None:
            raise AttributeError("R_values were never established in constructor")

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        var = 'Distance'

        ax1.set_title(var)
        ax1.set_xlabel('R')
        ax1.set_ylabel(var)

        set_r = list(set(self.rvals))
        set_r.sort(key=float, reverse=True)
        for i in particle_indexes:
            for n in ns:
                trans = []
                for r in set_r:
                    y = [AnalyzeNC(self.compute.traj.configurations[t].particles[i]).avg_distance_nth_particle(n)
                        for t in range(len(self.compute.traj.configurations))
                         if self.rvals[t] == r]
                    trans.append(np.average(y))
                ax1.plot(set_r, trans, label=str(i) + 'th nanoparticle ' + str(n) + ' in chain')
        plt.legend()
        plt.show()

    def orientation_v_timestep(self, particle_indexes):
        """

        :param particle_indexes: particle indexes to track orientations for
        :return: plot orientation on the surface of a sphere and coordinate angles theta and phi as a function of time
        """

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')

        var = 'Orientation'

        ax1.set_title(var)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        colors = ['g', 'brown', 'b', 'r', 'cyan', 'magenta'] * len(particle_indexes)
        for i in particle_indexes:
            ori = [AnalyzeCore(conf.particles[i].core).particle_1_spherical_coordinates() for conf in self.compute.traj.configurations]
            x = [np.sin(w[1]) * np.cos(w[0]) for w in ori]
            y = [np.sin(w[1]) * np.sin(w[0]) for w in ori]
            z = [np.cos(w[1]) for w in ori]
            ax1.plot(x, y, z, c=colors[i])
        plt.legend()
        plt.show()

        r = self.compute.traj.timesteps
        for i in particle_indexes:
            theta = [AnalyzeCore(conf.particles[i].core).particle_1_spherical_coordinates()[0] for conf in self.compute.traj.configurations]
            phi = [AnalyzeCore(conf.particles[i].core).particle_1_spherical_coordinates()[1] for conf in self.compute.traj.configurations]
            ax = plt.subplot(111, projection='polar')
            ax.plot(theta, r, label='theta ' + str(i))
            ax.plot(phi, r, label='phi' + str(i))
        ax.set_rmax(max(r))
        ax.grid(True)

        ax.set_title("Angles of orientation v timestep", va='bottom')
        plt.legend()
        plt.show()

    def heat_map(self, reference, particle_index=0 , dist_norm=True, angle=False):
        """
        creates a histogram measuring the relative differences in ligand structure between the trajectory particle and
        the reference particle

        :param reference: FunctionalizedParticle object to compare to. Must be the same core with same ligand.
        :param particle_index: the index of the particle from the trajectory to compare with
        :param dist_norm: When true v is calculated with adjustments for the differences in ligand length
        :param angle: when true the difference is given as angle in degrees, when false the difference is a value
        """

        v = []
        reference = cp.deepcopy(reference)
        reference.shift(np.multiply(-1, reference.core.position[0]))
        # reference.dump_xyz('reference')
        reference = reference.ligands
        for ind, con in enumerate(self.compute.traj.configurations):
            v_i = []
            part = cp.deepcopy(con.particles[particle_index])
            part.shift(- part.core.position[0])
            # part.dump_xyz(str(ind))
            for ind2, lig in enumerate(part.ligands):
                dif = AnalyzeLigand(lig).compare_pos(reference[ind2], dist_norm=dist_norm)
                v_i.append(dif)
            v.append(v_i)

        colors = ['g', 'brown', 'b', 'r', 'cyan', 'magenta'] * len(v)

        for i in range(len(v)):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            newdat = v[i]
            if angle:
                newdat = [np.degrees(np.arccos(-1 * ((v_sub_i**2 - 2) / 2))) for v_sub_i in v[i]]
                #print(newdat)
            binwidth = (max(newdat) - min(newdat)) / 20
            plt.hist(newdat, bins=np.arange(min(newdat), max(newdat) + binwidth, binwidth),
                     color=colors[i], alpha=0.5)
            ax.set_xlabel('V')
            if angle:
                ax.set_xlabel('Angle')
            ax.set_ylabel('counts')
            title = 'Chain Deformation Histogram '
            if self.rvals is not None:
                title += str(self.rvals[i])
            ax.set_title(title)
            plt.show()

    def three_d_map(self, reference, dist_norm=True, particle_index=0):
        """
        creates a 3d plot of relative ligand deformation as a function of the position of the first particle in the chain.

        :param reference: the reference FunctionalizedParticle object to compare to. Must be the same core with same ligand. Can be reinitialized from a simulation
        :param particle_index: the index of the particle from the trajectory to compare with
        :param dist_norm: When true v is calculated with adjustments for the differences in ligand length
        :param angle: when true the difference is given as angle in degrees, when false the difference is a value
        :return:
        """
        v = []
        reference2 = cp.deepcopy(reference)
        reference2.shift(np.multiply(-1, reference.core.position[0]))
        # reference.dump_xyz('reference')
        reference2 = reference2.ligands
        for ind, con in enumerate(self.compute.traj.configurations):
            v_i = []
            part = cp.deepcopy(con.particles[particle_index])
            part.shift(- part.core.position[0])
            # part.dump_xyz(str(ind))
            for ind2, lig in enumerate(part.ligands):
                dif = AnalyzeLigand(lig).compare_pos(reference2[ind2], dist_norm=dist_norm)
                v_i.append(dif)
            v.append(v_i)

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')

        var = 'Heat Map'

        ax1.set_title(var)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        for i in range(len(v)):
            for ind3 in range(len(reference.ligands)):
                x = reference.ligands[ind3].position[0][0]
                y = reference.ligands[ind3].position[0][1]
                z = reference.ligands[ind3].position[0][2]
                color = [v[0][ind3] / np.sqrt(2), 0, (np.sqrt(2) - v[0][ind3]) / np.sqrt(2)]
                ax1.scatter(x, y, z, c=color, s=1500)
            plt.legend()
            plt.show()

    def plot_density_grid(self, grid, rmax,dr,dtheta,dphi):
        """
        creates a 3d plot of the density of ligands as s function of spherical coordinates
        This plots the data created by the FunctionalizedConfiguration density_grid function so all the parameters should be the same

        :param grid: an array of data from the density_grid function of AnalyzeConfiguration
        :param rmax: max radius of the sphere
        :param dr: incremental r to divide the sphere into
        :param dtheta: incremental theta angle to divide the sphere into
        :param dphi: incremental phi angle to divide the sphere into
        """

        grid = np.sum(grid, axis=1)
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')

        var = 'Density'

        ax1.set_title(var)
        ax1.set_xlabel('r')
        ax1.set_ylabel('phi')
        ax1.set_zlabel('density')
        for i in range(grid.shape[1]):
            r = [dr * (.5 + w) for w in range(grid.shape[0])]
            phi = [dphi * (.5+ i) for w in range(grid.shape[0])]
            density = [grid[w][i] for w in range(grid.shape[0])]
            ax1.plot(r, phi, density)
        plt.legend()
        plt.show()

        ax1.set_title("Density", va='bottom')
        plt.legend()
        plt.show()

    def write_2d_mayavi_file(self, grid, rmax, dr,dtheta, dphi, filename='mayavi2.txt'):
        """
        writes a file used by a certain mayavi graphics script to do the density_grid using mayavi
        """
        grid = np.sum(grid, axis=1)
        f = open(filename,'w')
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                r = dr * (.5 + i)
                phi = dphi * (.5 + j)
                d = grid[i][j]
                f.write(str(r) + ' ' + str(phi) + ' ' + str(d) + '\n')

    def write_mayavi_file(self, reference, dist_norm=True, particle_index=0, filename='mayavi.txt'):
        """
        writes a file used by a certain mayavi graphics script to do 3d_map using mayavi
        """
        v = []
        reference2 = cp.deepcopy(reference)
        reference2.shift(np.multiply(-1, reference.core.position[0]))
        # reference.dump_xyz('reference')
        reference2 = reference2.ligands
        for ind, con in enumerate(self.compute.traj.configurations):
            v_i = []
            part = cp.deepcopy(con.particles[particle_index])
            part.shift(- part.core.position[0])
            # part.dump_xyz(str(ind))
            for ind2, lig in enumerate(part.ligands):
                dif = AnalyzeLigand(lig).compare_pos(reference2[ind2], dist_norm=dist_norm)
                v_i.append(dif)
            v.append(v_i)
        x = [reference.ligands[ind].position[0][0] for ind in range(len(v[0]))]
        y = [reference.ligands[ind].position[0][1] for ind in range(len(v[0]))]
        z = [reference.ligands[ind].position[0][2] for ind in range(len(v[0]))]

        f = open(filename, 'w')
        for i in range(len(x)):
            f.write(str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + ' ' + str(v[0][i]) + '\n')
