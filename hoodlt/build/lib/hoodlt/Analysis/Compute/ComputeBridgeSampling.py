"""
:module: Compute_BridgeSampling
:platform: Unix, Windows
:synopsis: computes free energies using MBAR

.. moduleauthor: Alex Travesset <trvsst@ameslab.gov> July 2023
.. history:
..
"""

import numpy as np
from scipy import optimize


class ComputeBridgeSampling:

    def __init__(self, num_windows, num_frames, bond_data, bond_types, k_spring, r_0, dict_optimize=None):
        """
        initializer for the class

        :param num_windows: number of windowa
        :param num_frames: number of configurations per window
        :param bond_data: numpy array with the spring displacements
        :param bond_types: parameterizes the type of bonds for each one of the third dimensions of the bond_data
        :param k_spring: values of spring for each window
        :param r_0: rest position for each window
        :param dict_optimize: parameters for optimization
        """

        self.num_files = num_windows
        self.num_frames = num_frames
        self.bond_data = bond_data
        self.dim_bonds = bond_types
        self.bond_k = k_spring
        self.bond_r0 = r_0
        if dict_optimize is None:
            self.dict_optimize = {'xtol': 1e-5}
        else:
            self.dict_optimize = dict_optimize

        self.solution = None

    def compute_solution(self, f_ini, overflow_val=500, u_fix_point=False):
        """
        returns the solution of the bridge sampling equation

        :param f_ini: initial values for the free energy
        :param overflow_val: the maximum value is :math:``\\exp(\\mbox{overflow_val})``
        :param u_fix_point: whether to use fixed_point or fsolve

        :return: ndarray the free energy computed by bridge sampling
        """

        if self.solution is None:
            self._solve_fp(f_ini, overflow_val=overflow_val, u_fix_point=u_fix_point)

        return self.solution

    def _uterm(self, i_term):
        """
        Returns the Harmonic energy calculated with the i-term

        :param i_term: index
        :return : ndarray where the argument represents the n-th configuration for the j-term
        """
        egy = np.zeros([self.num_files, self.num_frames])

        ini_val = 0
        for ind, ind_num in enumerate(self.dim_bonds):
            end_val = ini_val + ind_num
            mat_h = self.bond_k[i_term, ind]*(self.bond_data[:, :, ini_val:end_val]-self.bond_r0[i_term, ind])**2
            egy += 0.5*np.sum(mat_h, axis=2)
            ini_val = end_val

        return egy

    def _solve_fp(self, f_ini, overflow_val=500, u_fix_point=False):
        """

        :param f_ini: initial values for the free energy
        :param overflow_val: the maximum value is :math:``\\exp(\\mbox{overflow_val})``
        :param u_fix_point: whether to use fixed_point or fsolve

        :return: ndarray the free energy computed by bridge sampling
        """

        def fun_fp(f):
            return self.f_kernel(f, overflow_val=overflow_val)

        def fun_fs(f):
            return fun_fp(f)-f

        if u_fix_point:
            sol = optimize.fixed_point(fun_fp, f_ini, **self.dict_optimize)
        else:
            sol = optimize.fsolve(fun_fs, f_ini, **self.dict_optimize)

        self.solution = sol

    def f_kernel(self, f_arg, overflow_val=500):
        """
            Returns the kernel function of bridge sampling

            :param f_arg: initial values for the free energy
            :param overflow_val: the maximum value is :math:``\\exp(\\mbox{overflow_val})``

            :return: ndarray the free energy computed by bridge sampling
        """

        f_val = np.zeros([self.num_files, self.num_files, self.num_frames])
        for ind_l in range(self.num_files):
            f_val[ind_l, :, :] = f_arg[ind_l] - self._uterm(ind_l)

        f_max = np.max(f_val, axis=0)[np.newaxis, :, :]

        f_2 = np.where(f_max - f_val < overflow_val, f_val - f_max, -overflow_val)
        f_sum = f_max[0, :, :] + np.log(np.sum(np.exp(f_2), axis=0))
        for ind_i in range(self.num_files):
            f_val[ind_i, :, :] = f_sum[:, :] + self._uterm(ind_i)[:, :]

        f_min = np.min(f_val, axis=(1, 2))[:, np.newaxis, np.newaxis]
        f_1 = np.where(f_val - f_min < overflow_val, f_min - f_val, -overflow_val)

        val = f_min[:, 0, 0] - np.log(np.sum(np.exp(f_1), axis=(1, 2)))

        return val - val[-1]
