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
..                  - simplified code
..               Alex Travesset <trvsst@ameslab.gov? may 2023
..                  - Added Bridge sampling averages
..                  - Added function to compute potential energy and bias energy
..               Alex Travesset
..                  - solved the wham equation by bridge sampling
"""

import numpy
import numpy as np
from scipy import optimize
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Analysis.Collect.CollectLogData import CollectLogData
from hoodlt.Analysis.Collect.CollectHistData import CollectHistData
from hoodlt.Analysis.Collect.CollectForceFieldData import CollectForceFieldData
from hoodlt.Analysis.Collect.HistLogConsistency import consistency_therm_and_hist_files
from hoodlt.Analysis.Compute.ComputeBridgeSampling import ComputeBridgeSampling


class ComputeWhamPotential:
    """
    This class computes wham potentials, thermodynamic quantities, and potential of mean force from data outputted by
    2-body simulations
    """
    def __init__(self, file_names, ff_name, d_r=(0.2, 10), num_frames=None, end_val=None, qo=(-1, -1)):
        """

        :param file_names: names of the files outputted from the simulation, without extensions
        :param ff_name: name of the forcefield used to run the simulation
        :param d_r: width of the bins of the histograms that will be made
        :param num_frames: number of frames
        :param end_val: end value
        :param qo: return just two quantities (in the order stored) is useful to just return a subset of data
        """

        self.file_names = file_names
        self.ff_name = ff_name
        self.num_frames = num_frames
        self.end_val = end_val
        self.quant = 'potential_energy'
        self.delta_res = d_r[0]
        self.delta_quant = d_r[1]

        coll_hist = CollectHistData(self.file_names, num_frames=num_frames, end_val=end_val)
        coll_log = CollectLogData(self.file_names, self.ff_name, num_frames=num_frames, end_val=end_val, quant_only=qo)

        # ensure that the therm and hist files are taken at the exact same timesteps
        consistency_therm_and_hist_files(coll_log, coll_hist)

        bond_name = coll_hist.dict_bond_data['types']
        if len(bond_name) > 1:
            raise ValueError('PMF is calculated with a single spring. This file contains more than one. Stopping')
        self.bond_name = bond_name[0]
        self.lamda_all = coll_hist.lamda[:, 0]
        self.lamda = self.lamda_all[0]
        if not all([self.lamda_all[ind] == self.lamda for ind in range(coll_hist.num_files)]):
            print('Different lamda values accross the wham calculation')
        self.units = ForceFieldReader(ff_name).get_units()

        coll_ff = CollectForceFieldData(ff_name)

        # collect all the necessary simulation output data
        self.hist_data = coll_hist.bond_data[:, :, 0]  # can assume only one column
        self.simulation_temp = coll_hist.dict_bond_data['temperature']
        kbt = self.simulation_temp

        self.potential_energies = coll_log.get_quantity(self.quant)/kbt
        bk = coll_ff.get_bond_k(self.bond_name)/kbt

        # storing argument values
        self.bond_r0 = coll_hist.bond_r0[:, 0]
        self.bond_k = [bk*param for param in self.lamda_all]
        self.bond_data = coll_hist.bond_data
        self.rvals = self.bond_r0
        self.num_windows = len(self.rvals)

        # other data to store
        self.solved = None
        self.solution = None
        self.energies = None
        self.norm_rho_constant = None
        self.num_points = len(self.hist_data[0])
        self.hist, self.r_val = self._make_histograms(self.hist_data, self.delta_res)
        n_intervals = self.r_val.shape[0]-1
        # this is necessary to make sure indices are correctly assigned
        mat_ind = np.digitize(self.hist_data, self.r_val)-1
        mat_ind = np.where(mat_ind < 0, 0, mat_ind)
        self.hist_index = np.where(mat_ind > n_intervals-1, n_intervals-1, mat_ind)
        self.r = self.r_val[:-1] + 0.5 * np.diff(self.r_val)
        self.hist_2d = None
        self.u = None

        # necessary for intermediate calculation
        self.overflow_val = 600
        self._r = self.r[np.newaxis, :]
        self._k = np.array(self.bond_k)[:, np.newaxis]
        self._r0 = self.bond_r0[:, np.newaxis]

    def PMF_approx(self):
        """
        computes the PMF via integration

        :return: tuple (x, y) with points (x) and values of the pmf (y) at a given point
        """

        x = np.average(self.hist_data, axis=1)
        derivs = [self.bond_k[ind] * (dist - self.bond_r0[ind]) for ind, dist in enumerate(x)]

        g = [-np.trapz(derivs[x2:], x[x2:]) for x2 in range(len(x))]
        y = [g[-1] - g[x2] for x2 in range(len(x))]
        return x, y

    def solve(self, ini_f_val=None, d_optimize={'xtol': 1e-8}, over_flow_value=500, fp_solver=False):
        """
        Solves the wham calculation

        :param ini_f_val: initial configuration for the values of math:'f_i = \\Delta F_i'
        :param d_optimize: optimization parameters

        :return: None
        """
        if not self.solved:
            if ini_f_val is None:
                # estimate the values of f from approximate solution
                ini_conf = self._make_from_approx_pmf()
            else:
                # use provided initial value
                ini_conf = ini_f_val

            b_types = [1]
            b_k = np.array(self.bond_k)[:, np.newaxis]
            b_r0 = np.array(self.bond_r0)[:, np.newaxis]

            b_sampling = ComputeBridgeSampling(self.num_windows, self.num_frames, self.bond_data, b_types,
                                               b_k, b_r0, dict_optimize=d_optimize)

            self.solution = b_sampling.compute_solution(ini_conf, overflow_val=over_flow_value, u_fix_point=fp_solver)
            self.solved = True

    def solve_legacy(self, ini_f_val=None, xtol=1e-5, max_iterations=10000):
        """
        Does the wham calculation

        :param ini_f_val: initial configuration for the values of math:'f_i = \\Delta F_i'
        :param xtol: tolerance for solution
        :param max_iterations: maximum number of iterations
        :return: None
        """
        if not self.solved:
            # initial guess
            if ini_f_val is None:
                # estimate the values of f from approximate solution
                ini_conf = self._make_from_approx_pmf()
            else:
                # use provided initial value
                ini_conf = ini_f_val
            # tolerance
            self._fp_solve(ini_conf, xtol, max_iterations)
            self.solved = True

    def energy_potential_spring(self):
        """
        Computes the potential energy removing the spring, provides the spring energy

        :return: contributions for the potential energy and spring energyes
        :rtype: tuple
        """

        se = 0.5*self._k*(self.hist_data-self._r0)**2
        pe = self.potential_energies-se

        return pe, se

    @staticmethod
    def _make_histograms(data_file, delta_f):
        """
        makes histogram from the given data
k
        :param data_file: data array
        :param delta_f: bin width
        :return: function :math:'\\H_i(R)$'
        :rtype: numpy.ndarray list
        """

        d_min = np.amin(data_file)
        d_max = np.amax(data_file)
        htol = 0
        num = int((d_max - d_min) / delta_f)

        bin_delt = np.linspace(d_min, d_max, num, endpoint=True)

        h_fun = np.zeros([data_file.shape[0], bin_delt.shape[0] - 1])
        for ind in range(data_file.shape[0]):
            hist, bin_e = np.histogram(data_file[ind], bins=bin_delt)
            h_fun[ind] = hist

        if np.any(h_fun < htol):
            raise ValueError('Error in computing histogram')

        return h_fun, bin_delt

    def _make_histograms2d(self):
        """
        makes a 2d histogram from the given data

        """

        # subtract the potential energy from the data
        quant_data = self.energy_potential_spring()[0]

        d_min_u = np.amin(quant_data)
        d_max_u = np.amax(quant_data)
        num_u = int((d_max_u - d_min_u)/self.delta_quant)
        bin_u = np.linspace(d_min_u, d_max_u, num_u, endpoint=True)

        h_fun = np.zeros([self.hist.shape[0], self.r_val.shape[0] - 1, bin_u.shape[0] - 1])
        for ind in range(h_fun.shape[0]):
            hist, bin_edges_r, bin_edges_u = np.histogram2d(self.hist_data[ind], quant_data[ind],
                                                            bins=(self.r_val, bin_u))
            mat = np.sum(hist, axis=1)
            if not np.all(mat == self.hist[ind]):
                print('this error comes from histogram2d not compatible with histogram (1d)')
                print('check numpy histogram/histogram2 documentation')
                raise ValueError('Error in computing 2d histogram')
            h_fun[ind] = hist

        # this factor needs to be added
        self.hist_2d = h_fun
        self.u = bin_u[:-1] + 0.5 * np.diff(bin_u)

    def _fp_solve(self, ini_conf, xs, max_iterations):
        """
        solves the fixed point equation

        :param ini_conf: initial guess
        :param xs: tolerance
        :param max_iterations: maximum number of iterations
         """

        val = self._harm()
        min_val = np.min(val, axis=1)
        min_s = min_val[:, np.newaxis]
        f_eval = np.where(val-min_s > self.overflow_val, -self.overflow_val, -val+min_s)

        def fun_objective(x):
            rho = self.rho_f(x)
            iter_val = min_val - np.log(np.trapz(np.exp(f_eval)*rho, self.r, axis=1))
            return iter_val-iter_val[-1]

        self.solution = optimize.fixed_point(fun_objective, ini_conf, xtol=xs, maxiter=max_iterations)

    def _harm(self):
        """
        harmonic potential

        :return: Harmonic potential for each window at the center of the histogram values
        :rtype: numpy.ndarray
        """

        val = self._r - self._r0

        return 0.5*self._k*val*val

    def _make_from_approx_pmf(self):
        """
        approximate values for math:'f_i' based on the approximate pmf
        :return: numpy array
        """
        x, y = self.PMF_approx()
        pmf = np.array(y)[:, np.newaxis]
        val = self._harm()
        min_val = np.min(val, axis=1)
        min_s = min_val[:, np.newaxis]
        f_eval = np.where(val - min_s > self.overflow_val, -self.overflow_val, -val + min_s)
        f_approx = min_val - np.log(np.trapz(np.exp(f_eval-pmf), self.r, axis=1))
        return f_approx-f_approx[-1]

    def _intermediate_rho(self, x):
        """
        :math:'\\frac{1}{\\sum_{i=1}^N \\exp(-\beta (W_i(R)-f_i))} '

        :param x: :math:'\\Delta F_i'
        :return: value at x
        :rtype: tuple with the different parts of the calculation for later use
        """

        vals = x[:, np.newaxis] - self._harm()
        max_val = np.max(vals, axis=0)
        delta_val = np.where(max_val - vals > self.overflow_val, -self.overflow_val, vals - max_val)
        intg_a = 1 / np.sum(np.exp(delta_val), axis=0)
        min_val = np.min(max_val)
        delta_val = np.where(max_val - min_val > self.overflow_val, -self.overflow_val, min_val - max_val)

        return intg_a, max_val, delta_val, min_val

    def rho_f(self, x):
        """
        function :math:'\\rho(r)' as a function of parameter x

        :param x: :math:'\\Delta F_i'
        :return: value at x
        :rtype: numpy.ndarray
        """
        h_bias = np.sum(self.hist, axis=0)

        intg_a, max_val, delta_val, min_val = self._intermediate_rho(x)
        intg_c = h_bias*intg_a
        func = np.exp(delta_val)*intg_c
        f_val = np.trapz(func, self.r)
        return func/f_val

    def rho(self):
        """
        The actual function :math:'\\rho(r)'

        :return: :math:'\\rho(r)'
        :rtype: numpy.ndarray
        """

        self._check_solution()

        return self.rho_f(self.solution)

    def prob(self):
        """ probability that separation r between particles is in window i


        :return: math:'p_i(r)'
        :rtype: numpy.ndarray
        """
        self._check_solution()

        intg_a, m_v, d_v, m_m = self._intermediate_rho(self.solution)

        vals = self.solution[:, np.newaxis] - self._harm()
        delta_val = np.where(m_v - vals > self.overflow_val, -self.overflow_val, vals - m_v)

        func = intg_a * np.exp(delta_val)

        return func

    def pmf(self):
        """
        Computes the pmf

        :return: numpy.ndarray
        """
        self._check_solution()

        return -np.log(np.abs(self.rho()))

    def hist_simulation(self, ind):
        """
        Function :math:'\\rho_i(R)' as obtained from simulation

        :return: :math:'\\rho_i(R)'
        :rtype: numpy.ndarray
        """

        return self.hist[ind]/np.trapz(self.hist[ind], self.r)

    def hist_predicted(self, ind):
        """
        Function :math:'\\rho_i(R)' (computed biased distribution) to compare with simulation

        :return: :math:'\\rho_i(R)'
        :rtype: numpy.ndarray
        """

        self._check_solution()
        pot = self._harm()
        val = self.rho() * np.exp(-pot[ind]+self.solution[ind])
        return val/np.trapz(val, self.r)

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

    def add_energies(self):
        """
        Allows for calculation of internal energy and entropy
        Must solve the WHAM first

        :return: None
        """

        if not self.solved:
            raise ValueError("WHAM must be solved before this can be done")
        self._make_histograms2d()
        self.energies = True

    def rho_r_u(self):
        """
        Returns  function :math:'\\rho_i(R,u)'

        :return: numpy array
        """

        if not self.energies:
            raise ValueError('You need to add the function add_energies first')

        h_bias = np.sum(self.hist_2d, axis=0)
        intg_a, max_val, delta_val, min_val = self._intermediate_rho(self.solution)

        f_const = np.exp(delta_val) * intg_a
        func = f_const[:, np.newaxis]*h_bias
        f_val = np.trapz(np.sum(func, axis=1), self.r)

        return func / f_val

    def internal_energy(self):
        """
        Computes the internal energy :math:'U_i'

        :return: The internal energy
        :rtype: numpy.ndarray
        """

        rho_u = np.sum(self.u*self.rho_r_u(), axis=1)
        rho = self.rho()
        return rho_u/rho

    def entropies(self):
        """
        computes the entropy :math:'TS_i = F_i-U_i'

        :return: entropy
        :rtype: numpy.ndarray
        """
        return np.subtract(self.internal_energy(), self.free_energies())

    def bridge_sampling(self, num_windows, quant='potential_energy'):
        """
        Computes the average of a quantity using bridge sampling

        :param num_windows: use the Boltzmann distribution with bias potential i. (If i=-1 do not use any bias)
        :param quant: quantity to use
        :return: expectation value for the quantity
        :rtype: numpy.ndarray
        """

        self._check_solution()

        if quant == self.quant:
            quant_data = self.energy_potential_spring()[0]
        else:
            coll_log = CollectLogData(self.file_names, self.ff_name, num_frames=self.num_frames, end_val=self.end_val)
            quant_data = coll_log.get_quantity(quant)

        intg_a, m_v, d_v, m_m = self._intermediate_rho(self.solution)

        q_averages = np.zeros(num_windows.shape[0])

        val = self._harm()
        min_val = np.min(val, axis=1)
        min_s = min_val[:, np.newaxis]
        func = np.where(val - min_s > self.overflow_val, -self.overflow_val, -val + min_s)

        for ind, n_wind in enumerate(num_windows):
            if n_wind == -1:
                fac = 0.0
            else:
                fac = func[n_wind]

            prob = np.exp(d_v)*intg_a*np.exp(fac)/self.num_points
            prob_c = prob[self.hist_index]
            norm = np.sum(prob_c)

            q_averages[ind] = np.sum(quant_data*prob_c)/norm

        return q_averages
