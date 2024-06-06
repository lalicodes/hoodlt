"""
:module: Compute_Lattice_Potential
:platform: Unix, Windows
:synopsis: computes potential of a superlattice  fcc, bcc, planar hexagonal, planar square, 1D chain.

.. moduleauthor: Xun Zha <xzha@iastate.edu> November 2018
.. history:
..                Xun Zha <xzha@iastate.edu> June 2021
..                  - fixed some errors
..                Xun Zha <xzha@iastate.edu> July 2021
..                  - bug fix
..                Alex Travesset <trvsst@ameslab.gov> June 2022
..                  - redid the class to incorporate all modifications from hoomd v3
"""

import numpy as np
from scipy import interpolate
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Data.Modelconfigurations.Saver import load_config
from hoodlt.Analysis.Collect.CollectLogData import CollectLogData
from hoodlt.Analysis.Collect.CollectHistData import CollectHistData
from hoodlt.Analysis.Collect.CollectForceFieldData import CollectForceFieldData
from hoodlt.Analysis.Collect.HistLogConsistency import consistency_therm_and_hist_files
from hoodlt.Analysis.Compute.ComputeBridgeSampling import ComputeBridgeSampling
import hoodlt.Utils.LatticeNeighbors as Ln


class ComputeFreeEnergyLattice:

    def __init__(self, file_names, ff_name,  u_d=1, num_frames=None, end_val=None, dict_optimize=None):
        """
        initializer for the class

        :param file_names: list of filenames
        :param ff_name: name of the forcefield
        :param u_d: reference distance (in construction units) suitable to express dimensionless pressure
        :param num_frames: number of frames
        :param end_val: end value
        :param dict_optimize: parameters for optimization
        """

        self.file_names = file_names
        self.num_files = len(self.file_names)
        self.ff_name = ff_name
        self.num_frames = num_frames
        self.end_val = end_val

        # obtain the data
        coll_log = CollectLogData(self.file_names, self.ff_name, num_frames=num_frames, end_val=end_val)
        self.num_frames = coll_log.num_frames
        self.end_val = end_val
        coll_hist = CollectHistData(self.file_names, num_frames=self.num_frames, end_val=self.end_val)
        coll_ff = CollectForceFieldData(ff_name)

        # ensure that the therm and hist files are taken at the exact same timesteps
        consistency_therm_and_hist_files(coll_log, coll_hist)

        # get the units
        self.units = ForceFieldReader(ff_name).get_units()
        self.unit_distance = self.units.length_construction_to_simulation*u_d

        # lattice object and number of (nano)particles in the lattice
        self.lattices = self.num_files*[0]
        self.vol_pp = np.zeros(self.num_files)
        for ind, names in enumerate(self.file_names):
            conf = load_config(names)
            self.lat = conf.lat
            self.lattices[ind] = self.lat
            # volume per particle
            self.vol_pp[ind] = self.lat.vol_unit_cell() / (self.unit_distance ** 3 * np.sum(self.lat.typ))

        self.num_lattice_points = self.lat.num_pnts()

        # get the temperature
        self.simulation_temp = coll_hist.dict_bond_data['temperature']
        kbt = self.simulation_temp

        # get bond data
        self.bond_name = coll_hist.dict_bond_data['types']
        self.lamda = coll_hist.lamda
        self.bond_r0 = coll_hist.bond_r0
        self.bond_k = np.array([[coll_ff.get_bond_k(bnd)*self.lamda[ind2, ind1]/kbt for ind1, bnd
                                 in enumerate(self.bond_name)] for ind2 in range(self.num_files)])
        self.dim_bonds = coll_hist.dim_bonds
        if isinstance(self.dim_bonds, np.ndarray) or isinstance(self.dim_bonds, list):
            pass
        else:
            self.dim_bonds = [1]
        self.bond_data = coll_hist.bond_data

        # store pressure and internal energy
        self.instant_pressure = self.unit_distance**3*coll_log.get_quantity('pressure')/kbt
        self.instant_int_energy = coll_log.get_quantity('potential_energy')/kbt

        # pressure
        self.pressure = np.average(self.instant_pressure, axis=1)
        # average internal energy per particle# lamda figure of merit
        self.int_energy_pp = np.average(self.instant_int_energy, axis=1)/np.sum(self.num_lattice_points)

        # lamda figure of merit
        f_l = 0.1
        f_n = f_l*self.unit_distance**3
        self.lamda_optimal = -self.pressure*self.bond_r0[:, 0]*self.lamda[:, 0]/(f_n*self.bond_k[:, 0])

        # bridge sampling
        self.bridge_sampling_solution = None
        self.dict_optimize = dict_optimize
        if dict_optimize is None:
            self.dict_optimize = {'xtol': 1e-8}
        ndim1, ndim2 = self.bond_r0.shape
        self._b = np.zeros([ndim1, ndim2, 1, 1])
        self._k = np.zeros_like(self._b)
        self._b[:, :, 0, 0] = self.bond_r0[:, :]
        self._k[:, :, 0, 0] = self.bond_k[:, :]

    def compute_avg_distance(self):
        """
        Provides the average distance for each bond

        :return : ndarray
        """

        avgs = np.zeros([self.num_files, len(self.bond_name)])
        ini_val = 0
        for ind, ind_num in enumerate(self.dim_bonds):
            end_val = ini_val + ind_num
            avgs[:, ind] = np.average(self.bond_data[:, :, ini_val:end_val], axis=(2, 1))
            ini_val = end_val

        return avgs

    def compute_hist_distances(self, nbins=10):
        """

        :param nbins: number of bins to use
        Provides the histogram for distances within the lattice

        :return: list of histograms
        """
        lst_all = []
        for ind in range(self.num_files):
            lst = []
            ini_val = 0
            for ind_num in self.dim_bonds:
                end_val = ini_val + ind_num
                lst.append(np.histogram(self.bond_data[ind, :, ini_val:end_val].flatten(), bins=nbins, density=True))
                ini_val = end_val
            lst_all.append(lst)

        return lst_all

    def compute_pressure_correction(self):
        """
        Computes the pressure

        :return: ndarray, the pressure including the correction from springs
        """

        c_sum = np.zeros([self.num_files])

        for i_term in range(self.num_files):
            ini_val = 0
            for ind, ind_num in enumerate(self.dim_bonds):
                end_val = ini_val + ind_num
                arg_r = self.bond_data[i_term, :, ini_val:end_val] - self.bond_r0[i_term, ind]
                mat_h = self.bond_k[i_term, ind] * self.bond_r0[i_term, ind]*np.sum(arg_r, axis=1)
                c_sum[i_term] += np.average(mat_h, axis=0)
                ini_val = end_val

        fac = 3.0*np.sum(self.num_lattice_points)*self.vol_pp[:]

        return c_sum/fac

    def compute_free_energy(self, with_correction_term=True):
        """
        Computes the free energy
        :param with_correction_term: Include the correction term

        :return: ndarray, the free energy computed by integrating the pressure
        """

        press = np.zeros_like(self.pressure)
        press[:] = self.pressure[:]
        if with_correction_term:
            press += self.compute_pressure_correction()

        pchip_integrand = interpolate.PchipInterpolator(self.vol_pp, press)
        v0 = self.vol_pp[-1]
        free_energy = [pchip_integrand.integrate(vol, v0) for vol in self.vol_pp]

        return free_energy

    def compute_free_energy_bridge_sampling(self, f_ini=None, over_flow_value=500, fp_solver=False, eval_kernel=False):
        """
        :param f_ini: initial values for the free energy
        :param over_flow_value: overflow value for exponential
        :param fp_solver: if True use fixed_point solver, if False use fsolve
        :param eval_kernel: evaluate the bridge sampling kernel with the initial condition f_ini. Does not solve

        :return: ndarray, the free energy computed by bridge sampling
        """

        b_sampling = ComputeBridgeSampling(self.num_files, self.num_frames, self.bond_data, self.dim_bonds,
                                               self.bond_k, self.bond_r0, dict_optimize=self.dict_optimize)
        of = over_flow_value
        fp = fp_solver
        if f_ini is None:
            i_c = np.array(self.compute_free_energy()) * np.sum(self.num_lattice_points)
        else:
            i_c = f_ini

        if eval_kernel:
            result = b_sampling.f_kernel(f_ini)
        else:
            if self.bridge_sampling_solution:
                result = b_sampling.compute_solution(i_c, overflow_val=of, u_fix_point=fp)
                self.bridge_sampling_solution = result
            else:
                result = self.bridge_sampling_solution

        return result/np.sum(self.num_lattice_points)

    def compute_free_energy_pairwise(self, e_function, degree=2):
        """
        Computes the free energy assuming that it is given as a sum of pairwise potentials

        :param e_function: pair potential
        :param degree: 1 for single component, 2 for binary, etc..
        :param degree: include 1=nearest neighbor, 2= next to nearest neighbor, etc..

        :return: ndarray, prediction of free energy for the given lattices
        """

        e_vals = np.zeros(len(self.lattices))

        for ind_lat, lat in enumerate(self.lattices):
            lat_neighb = Ln.LatNeighbor(lat)
            num_base = lat_neighb.base
            num_type = lat_neighb.typ
            num_pnts = len(num_base)
            egy = 0.0
            for ind1 in range(num_pnts):
                type_1 = num_type[ind1]
                for deg in range(1, degree):
                    nn_p = lat_neighb.neighbor(ind1, deg)
                    dist = lat_neighb.u_dist[num_base[ind1]][deg]
                    num_neighbors = len(nn_p)
                    type_2 = num_type[nn_p[0]]
                    egy += num_neighbors*e_function(type_1, type_2, dist)
            e_vals[ind_lat] = 0.5*egy/num_pnts

        return e_vals

    @staticmethod
    def _fit_point(x, y, fit_x):
        """PChip fit of y at given point of x

        :param x:
        :param y:
        :param fit_x:
        :return: fit_y
        """

        return interpolate.pchip_interpolate(x, y, fit_x)
