""":module: ComputeWhamPotential
   :platform: Unix, Windows
   :synopsis: computes the potential

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Moved old code (written by Alex Travesset) to fit new analysis scheme
..                Matthew Pham <mpham@iastate.edu>, May 2020
..                  - Added option to use different bond scaling for each rval/window
..                Alex Travesset <trvsst@ameslab.gov>, February 2022
..                  - Fixed documentation issues
..                  - Added the probability of p_i
..                  - Refined the calculation for large r tails and energy
..               Alex Travesset <trvsst@ameslab.gov>, May 2022
..                  - made it compatible with the new HOODLT
..
"""

import numpy as np
from scipy import trapz, optimize
from hoodlt.Analysis.Collect.CollectLogData import CollectLogData
from hoodlt.Analysis.Collect.CollectHistData import CollectHistData
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Analysis.Collect.CollectForceFieldData import CollectForceFieldData


class ComputeWhamPotential:
    """
    This class computes wham potentials, thermodynamic quantities, and potential of mean force from data outputted by
    2-body simulations
    """
    def __init__(self, file_names, ff_name, delta_res=0.2, num_frames=None, end_val=None, quant='potential_energy'):
        """

        :param file_names: names of the files outputted from the simulation, without extensions
        :param ff_name: name of the forcefield used to run the simulation
        :param delta_res: width of the bins of the histograms that will be made
        :param num_frames: number of frames
        :param end_val: end value
        :param quant: quantity to plot with the histogram
        """

        self.file_names = file_names
        self.ff_name = ff_name
        self.num_frames = num_frames
        self.end_val = end_val
        self.quant = quant

        coll_hist = CollectHistData(self.file_names, num_frames=num_frames, end_val=end_val)
        coll_log = CollectLogData(self.file_names, self.ff_name, num_frames=num_frames, end_val=end_val)

        bond_name = coll_hist.dict_bond_data['types']
        if len(bond_name) > 1:
            raise ValueError('PMF is calculated with a single spring. This file contains more than one. Stopping')
        self.bond_name = bond_name[0]
        self.lamda_all = coll_hist.lamda[:, 0]
        self.lamda = self.lamda_all[0]
        if not all([self.lamda_all[ind] == self.lamda for ind in range(coll_hist.num_files)]):
            raise ValueError('all values of lamda have to be the same for wham calculation')
        self.units = ForceFieldReader(ff_name).get_units()

        coll_ff = CollectForceFieldData(ff_name)

        # collect all the necessary simulation output data
        self.hist_data = coll_hist.bond_data[:, :, 0]  # can assume only one column
        self.simulation_temp = coll_hist.dict_bond_data['temperature']
        kbt = self.simulation_temp

        self.potential_energies = coll_log.get_quantity(quant)/kbt
        bk = coll_ff.get_bond_k(self.bond_name)/kbt

        # storing argument values
        self.bond_r0 = coll_hist.bond_r0[:, 0]
        self.bond_k = [bk*self.lamda]*len(self.bond_r0)
        self.rvals = self.bond_r0
        self.delta_res = delta_res
        self.num_windows = len(self.rvals)

        # other data to store
        self.solved = False
        self.solution = None
        self.wham2 = None
        self.energies = None
        self.norm_rho_constant = None
        self.num_points = len(self.hist_data[0])
        self.hist, r_val = self._make_histograms(self.hist_data, delta_res)
        self.r = r_val[:-1] + 0.5 * np.diff(r_val)

    def _make_histograms(self, data_file, delta_f):
        """
        computes the function :math: 'H_i(R)'

        :param data_file: data array
        :param delta_f: bin width
        :return: function :math:'$H_i(R)$'
        :rtype: numpy.ndarray list
        """

        d_min = np.amin(data_file)
        d_max = np.amax(data_file)
        # htol = 20
        htol = 0
        num = int((d_max - d_min) / delta_f)

        bin_delt = np.linspace(d_min, d_max, num, endpoint=True)

        h_fun = np.zeros([data_file.shape[0], bin_delt.shape[0] - 1])

        for ind in range(h_fun.shape[0]):
            hist, bin_edges = np.histogram(data_file[ind], bins=bin_delt)
            h_fun[ind] = hist

        for i in range(len(h_fun[ind])):
            if h_fun[ind][i] < htol:
                h_fun[ind][i] = 0.0

        # print(bin_edges)
        return h_fun, bin_edges

    def _w_int(self, xval, ind):
        """
        Helper function, convenient for direct iterations and checking solutions

        :param xval: x input
        :param ind: window number
        :return: integrand value
        :rtype:    numpy.ndarray
        """

        return trapz(self._int_rho(xval) * self._harm(ind), self.r)

    def _harm(self, ind):
        """
        harmonic potential

        :param ind: window number
        :return: Harmonic potential at the center of the histogram values
        :rtype: numpy.ndarray
        """

        val = (self.r - self.bond_r0[ind]) * (self.r - self.bond_r0[ind])

        return np.exp(-0.5 * self.bond_k[ind] * val)

    def PMF_approx(self):
        """
        computes the approximate PMF via integration if the wham cannot be solved

        :return: a tuple x, y
        """

        x = np.average(self.hist_data, axis=1)
        derivs = [self.bond_k[i] * (x[i] - self.bond_r0[i]) for i in range(len(x))]

        g = [trapz(derivs[:x2], x[:x2]) for x2 in range(len(x))]
        y = [g[0] - g[x2] for x2 in range(len(x))]
        return x, y

    def _fp_solve(self, ini_conf, xs=1e-8):
        """
        Computes the quantity :math: '

        :param ini_conf: initial guess
        :param xs: tolerance
        """

        self.solution = np.zeros([self.num_windows])
        val_fun = np.zeros([self.num_windows])

        def loc_fun(x):
            for ind in range(self.num_windows):
                val_fun[ind] = 1.0 / self._w_int(x, ind)
            return val_fun

        self.solution = optimize.fixed_point(loc_fun, ini_conf, xtol=xs, maxiter=10000)

    def _int_rho(self, x):
        """
        function :math:'\\rho(r)' as a function of parameter x

        :param x: :math:'\\exp(\\beta \\Delta F_i)'
        :return: value at x
        :rtype: numpy.ndarray
        """

        norm = np.zeros_like(self.r)

        for ind in range(self.num_windows):
            norm += self.num_points * self._harm(ind) * x[ind]

        val = np.sum(self.hist, axis=0) / norm

        self.norm_rho_constant = trapz(val, self.r)

        return val / self.norm_rho_constant

    def rho(self):
        """
        The actual function :math:'\\rho(r)'

        :return: :math:'\\rho(r)'
        :rtype: numpy.ndarray
        """

        self._check_solution()

        return self._int_rho(self.solution)

    def pmf(self):
        """
        Computes the pmf

        :return: numpy.ndarray
        """
        self._check_solution()

        return -np.log(np.abs(self.rho()))

    def h_original(self, ind):
        """
        Function :math:'\\rho_i(R)' as obtained from simulation

        :return: :math:'\\rho_i(R)'
        :rtype: numpy.ndarray
        """

        int_norm = np.zeros(self.num_windows)
        for ind_p in range(self.num_windows):
            int_norm[ind_p] = trapz(self.hist[ind_p], self.r)

        return self.hist[ind] / int_norm[ind]

    def h_comp(self, ind):
        """
        Function :math:'\\rho_i(R)' (computed biased distribution) to compare with simulation

        :return: :math:'\\rho_i(R)'
        :rtype: numpy.ndarray
        """

        self._check_solution()

        return self._harm(ind) * self.rho() * self.solution[ind]

    def prob(self, ind):
        """ probability that separation r between particles is in window i


        :return: math:'p_i(r)'
        :rtype: numpy.ndarray
        """

        return self.pre_prob()*self._harm(ind) * self.solution[ind]

    def pre_prob(self):
        """

        :return:  math:'e^{\\beta (W_i(R)-F_i)} p_i'
        :rtype:  numpy.ndarray
        """

        self._check_solution()

        norm = np.zeros_like(self.r)

        for ind in range(self.num_windows):
            norm += self.num_points * self._harm(ind) * self.solution[ind]

        return self.num_points / norm

    def _check_solution(self):
        """
        helper function checking whether the wham has already been solved

        :return: None
        """

        if type(self.solution).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
        return 0

    def free_energies(self):
        """
        Returns the free energy normalized to zero

        :return: :math:'F_i'
        :rtype: numpy.ndarray
        """

        f_i = self.pmf() - self.pmf()[-1]

        return f_i

    def solve(self, ini_f_val=None):
        """
        Does the wham calculation, after the wham is solved, it can be plotted and other thermodynamic quantities can
        be calculated

        :param ini_f_val: initial configuration for the values of math:'f_i = e^{\\beta \\Delta F_i)} '
        :return: None
        """
        if not self.solved:
            # initial guess
            if ini_f_val is None:
                # ini_conf = 0.8*np.ones(self.num_windows) we used this value in the past
                # estimate the values of f from approximate solution
                ini_conf = self._make_from_approx_pmf()
                print(ini_conf)
            else:
                # use provided initial value
                ini_conf = ini_f_val
            # tolerance
            xtol = 1e-8
            self._fp_solve(ini_conf, xtol)
            self.solved = True

    def _make_from_approx_pmf(self):
        """
        Provides approximate values for math:'f_i' based on the approximate pmf
        :return: numpy array
        """
        x, y = self.PMF_approx()
        rv = self.r
        fnc = np.zeros_like(rv)
        for ind in range(fnc.shape[0]):
            indx = np.abs(x - rv[ind]).argmin()
            fnc[ind] = y[indx]
        rho_approx = np.exp(-fnc)
        cons = np.trapz(rho_approx, rv)
        f_approx = np.zeros(self.num_windows)
        for ind in range(self.num_windows):
            f_approx[ind] = cons / np.trapz(self._harm(ind) * rho_approx, rv)
        return f_approx

    def add_energies(self, delta_res_energy=10):
        """
        Allows for calculation of internal energy and entropy
        Must solve the WHAM first
        Only works if internal energy was logged at the same rate that the distances between the particles were sampled

        :param delta_res_energy: bin width for energy histograms
        :return: None
        """

        if not self.solved:
            raise ValueError("WHAM must be solved before this can be done")
        wham2 = self.make_wham_energy(self.potential_energies, delta_res_energy)
        self.wham2 = wham2
        self.energies = True

    def make_wham_energy(self, data_file, delta_f):
        """
        computes the function :math: '$H_i(R)$'

        :param data_file: data array of the energies
        :param delta_f: bin width for the energies
        :return: 2d WhamPotential object
        :rtype: WhamPotential2d
        """

        wham_energy = ComputeWhamPotential2D(self, data_file, delta_f)

        return wham_energy

    def internal_energy(self):
        """
        Computes the internal energy :math:'U_i'

        :return: The internal energy
        :rtype: numpy.ndarray
        """
        if self.wham2 is None:
            raise ValueError("Must call add2d before using this function, internal_energy")

        return self.wham2.internal_energy()

    def entropies(self):
        """
        computes the entropy :math:'TS_i = F_i-U_i'

        :return: entropy
        :rtype: numpy.ndarray
        """
        return np.subtract(self.internal_energy(), self.free_energies())


class ComputeWhamPotential2D(ComputeWhamPotential):
    """
    wham potential with a 2d histogram for computing the internal energy
    """

    def __init__(self, orig_wham, energy_data_file, delta_f_energy):
        """
        Constructor

        :param orig_wham: wham object with :math:'\\rho(r)' computed
        :param energy_data_file: a numpy array  [num_of_windows, number of data points per window]
        :param delta_f_energy: the width of the bins
        """
        super(ComputeWhamPotential2D, self).__init__(orig_wham.file_names, orig_wham.ff_name, orig_wham.delta_res,
                                                     orig_wham.num_frames, orig_wham.end_val, orig_wham.quant)

        self.solution = orig_wham.solution

        if energy_data_file.shape[0] != self.num_windows:
            raise ValueError("energy data must have same number of windows")

        if energy_data_file.shape[1] != self.num_points:
            raise ValueError("energy data must have same number of points in each window")

        for ind in range(energy_data_file.shape[0]):
            for ind2, u in enumerate(energy_data_file[ind]):
                energy_data_file[ind][ind2] = u - self.harmonic_energy(ind, self.hist_data[ind][ind2])
        d_min = np.amin(self.hist_data)
        d_max = np.amax(self.hist_data)
        num = int((d_max - d_min) / self.delta_res)
        bin_delt = np.linspace(d_min, d_max, num, endpoint=True)

        d_min_u = np.amin(energy_data_file)
        d_max_u = np.amax(energy_data_file)
        num_u = int((d_max_u - d_min_u) / delta_f_energy)
        bin_delt_u = np.linspace(d_min_u, d_max_u, num_u, endpoint=True)

        h_fun = np.zeros([energy_data_file.shape[0], bin_delt.shape[0] - 1, bin_delt_u.shape[0] - 1])

        for ind in range(h_fun.shape[0]):
            hist, bin_edges_r, bin_edges_u = np.histogram2d(self.hist_data[ind], energy_data_file[ind],
                                                            bins=(bin_delt, bin_delt_u))
            h_fun[ind] = hist

        self.r = bin_edges_r[:-1] + 0.5 * np.diff(bin_edges_r)

        self.u = bin_edges_u[:-1] + 0.5 * np.diff(bin_edges_u)

        self.hist = h_fun

        self.o_wham = orig_wham

    def _int_rho(self, x):
        """
        function :math:'$\\rho$'

        :param x: :math:'$\\exp(\\beta \\Delta F_i)$'
        :return: value at x
        :rtype: numpy array                            elf.rho()
        """

        norm = np.zeros_like(self.r)

        for ind in range(self.num_windows):
            norm += self.num_points * self._harm(ind) * x[ind]

        norm2 = np.zeros((self.hist.shape[1], self.hist.shape[2]))

        for x in range(norm2.shape[0]):
            for i in range(norm2.shape[1]):
                norm2[x][i] = norm[x]

        norm = norm2

        val = np.sum(self.hist, axis=0) / norm

        val_const = trapz(trapz(val, self.u), self.r)
        # normalize the result

        val = val / val_const

        # print("should be 1")
        # print(trapz(trapz(val, self.u), self.r))
        return val

    def internal_energy(self):
        """
        computes :math:'U(r)' the internal energy at parameter r

        :return: internal energy at parameter r
        :rtype: numpy.ndarray
        """

        h_fun_u = np.zeros([self.r.shape[0], self.u.shape[0]])

        for ind in range(self.num_windows):
            h_fun_u += self.hist[ind]

        h_fun_u = np.transpose(np.transpose(h_fun_u)*self.o_wham.pre_prob())
        norm = trapz(trapz(h_fun_u, self.u, axis=1), self.r)
        uval = trapz(h_fun_u*self.u, self.u, axis=1)

        return uval/(norm*self.o_wham.rho())

    def internal_energy_old(self):
        """
        Returns the internal energy of the system  (old calculation)

        :return: the internal energy as a function of r
        """
        rho = self.rho()
        u_r = np.zeros(rho.shape[0])

        for ind in range(u_r.shape[0]):
            u_r[ind] = trapz(rho[ind] * self.u, self.u)

        return u_r / self.o_wham.rho()

    def normal_rho(self):
        """

        :return: rho should be proportional to the original rho; off by a factor of the ratio of integral normalizations
        """

        x = self.rho()
        rho_r = trapz(x, self.u)
        # print("should be 1.0")
        # print(trapz(self.o_wham.rho(), self.r))
        return rho_r

    def harmonic_energy(self, win_num, r):
        """
        Average energy of the harmonic bond in the given window

        :param win_num: window number
        :param r: r value
        :return: harmonic bias enegy
        """

        val = (r - self.bond_r0[win_num]) * (r - self.bond_r0[win_num])
        return 0.5 * self.bond_k[win_num] * val
