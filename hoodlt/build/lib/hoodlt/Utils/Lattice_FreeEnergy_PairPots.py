"""
:module: BinaryLattice_FreeEnergy_PairPots
:platform: Unix, Windows
:synopsis:

.. moduleauthor: Nathan Horst <nhorst@iastate.edu> March 2020
.. history::
"""

from __future__ import division
import copy as cp
import numpy as np
from hoodlt.Utils.Lattice_UnitCell import unit_cell as Uc
import numpy.linalg as la
from scipy import stats
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from hoodlt.Utils.LatticeNeighbors import LatNeighbor as Ln
import gsd.hoomd
from hoodlt.Data.Modelconfigurations.Saver import load_config
import hoodlt.Lattices.reciprocal_vectors as rv
import hoodlt.Groups.Quaternion as Quaternion


class BinaryLatticeFreeEnergyPairPots(object):

    def __init__(self, lat, pair_pots):
        """
        :param lat: Binary Lattice Object
        :param pair_pots: list of pair potential filenames -- in the order [A-B, A-A, B-B]
        """

        if len(pair_pots) != 3:
            raise ValueError('Input A-B, A-A, B-B interaction potential')

        r_val_ab, pair_pot_ab = np.loadtxt(pair_pots[0], unpack=True)
        r_val_aa, pair_pot_aa = np.loadtxt(pair_pots[1], unpack=True)
        r_val_bb, pair_pot_bb = np.loadtxt(pair_pots[2], unpack=True)

        pair_pot_ab -= pair_pot_ab[-1]
        pair_pot_aa -= pair_pot_aa[-1]
        pair_pot_bb -= pair_pot_bb[-1]

        self.ri = [r_val_ab, r_val_aa, r_val_ab, r_val_bb]
        self.pi = [pair_pot_ab, pair_pot_aa, pair_pot_ab, pair_pot_bb]
        self.lat = lat
        #self.unit_cell = Uc(self.lat)
        self.lnn = Ln(self.lat)
        self.calculated = False

    def calc_free_energy(self, a_nn_range, res=50):
        """
        :param a_nn_range: range of a_nn values for free energy calculation (list)
        :return: free energy for a_nn
        """

        self.calculated = True
        self.pp = np.zeros(res)
        self.a_nn_range = np.array(a_nn_range)
        if self.a_nn_range[-1] < self.a_nn_range[0]:
            raise ValueError('Reverse list order...')

        self.a_ratio = np.linspace(self.a_nn_range[0]/self.a_nn_range[-1], 1, res)

        ind_a = np.where(self.lnn.typ == 0)[0][0]
        ind_b = np.where(self.lnn.typ == 1)[0][0]
        # search degrees for A-B, A-A, B-A, B-B
        if np.all(self.lnn.typ[self.lnn.neighbor(ind_a, 1)]):  # first nearest neighbors are of type B
            degree_a = np.array([1, 2])
        else:  # first nearest neighbors are of type A
            degree_a = np.array([2, 1])
        if np.all(self.lnn.typ[self.lnn.neighbor(ind_b, 1)]):  # first nearest neighbors are of type B
            degree_b = np.array([2, 1])
        else:  # first nearest neighbors are of type A
            degree_b = np.array([1, 2])
        indices = np.array([ind_a, ind_a, ind_b, ind_b])
        degrees = np.concatenate((degree_a, degree_b))
        num_typs = np.array([self.lat.typ[0], self.lat.typ[0], self.lat.typ[1], self.lat.typ[1]])
        self.pplist = []
        for i in range(4):
            print(self.lnn.u_dist[indices[i]][degrees[i]], indices[i], degrees[i], '\n')
            a_val_i = self.a_ratio*self.lnn.u_dist[indices[i]][degrees[i]]
            #print(a_val_i, '\n')
            pp_i = self.pi[i]*len(self.lnn.neighbor(indices[i], degrees[i]))*num_typs[i]
            self.pp += np.interp(a_val_i, self.ri[i], pp_i)
            self.pplist.append(np.interp(a_val_i, self.ri[i], pp_i))
            #print(np.interp(a_val_i, self.ri[i], pp_i))
        self.pp /= np.sum(num_typs)

        return self.pp


    def plot_free_energy(self,file_n, title='Free Energy from Pair Potentials', xlabel='a$_{nn}$', ylabel='Free Energy'):
        """
        :param file_n: file name of plot
        :param title: title of plot
        :param xlabel: label for x axis
        :param ylabel: label for y axis
        :return:
        """

        if not self.calculated:
            raise ValueError('Must calculate Free Energy before plotting.')
        rvals = np.linspace(self.a_nn_range[0], self.a_nn_range[-1], len(self.a_ratio))
        np.savetxt(file_n+'.pmf', np.vstack((rvals, self.pp)).T)

        plt.rc('text', usetex=True)
        plt.rc('font', family = 'serif', serif ='cm10')

        plt.title(title, fontsize=18)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)

        ppf = PdfPages(file_n+'_free_energy.pdf')
        #ppf.savefig()


        print(rvals[np.argmin(self.pp)], 'Rmin')
        print(np.amin(self.pp), 'FE min')
        print(rvals[np.argmin(self.pplist[0])] - rvals[np.argmin(self.pp)], 'deltar')
        print(np.amin(self.pp) - np.amin(self.pplist[0]), 'deltaE')
        print(self.pplist[1][np.argmin(self.pplist[0])], 'repuls at min')

        plt.plot(rvals, self.pp, label= 'PP')
        plt.plot(rvals, self.pplist[0], label = 'A-B')
        plt.plot(rvals, self.pplist[1], label = 'A-A')
        plt.plot(rvals, self.pplist[2], label = 'B-A')
        plt.plot(rvals, self.pplist[3], label = 'B-B')

        plt.legend()
        # save figure to pdf
        ppf.savefig()
        # close PdfPages
        ppf.close()
        plt.close()

