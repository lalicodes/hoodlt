"""
:module: D_matrix_Mixt_D
:platform: Unix, Windows
:synopsis: Defines the class for the calculation of the D-matrix

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December2015
"""

import numpy as np

import hoodlt.D_matrix.Dmatrix_Mixt_distances as Dm
from numpy.fft import fftn


class DMatrix(object):
    """ Defines the class that actually computes the Lattice sums and the D-matrix

    """

    def __init__(self, lat, pot):
        """The constructor

        :param lat: lattice object
        :param pot: potential object
        """

        self.pot = pot
        self.lat = lat
        self.d_d = Dm.MinConv(lat, pot.cut)

        # number of basis
        n_bas = np.sum(lat.typ)
        # number of species
        n_spc = len(lat.typ)

        self.n_spc = n_spc
        self.num_species = float(n_bas)
        self.ncoord = np.zeros(n_spc, dtype=int)
        self.lcoord = np.zeros(n_spc, dtype=int)

        # define the internal coordinates
        bfac = 0
        for dims in range(n_spc):
            self.ncoord[dims] = bfac
            bfac += lat.typ[dims]
            self.lcoord[dims] = bfac

        # r-vec is a matrix stores that stores all lattice positions
        self.r_vec = np.zeros([3, n_bas, n_bas, lat.l[0], lat.l[1], lat.l[2]])

        # D-matrix in real and k-space
        self.d_real = np.zeros([3*n_bas, 3*n_bas, lat.l[0], lat.l[1], lat.l[2]])
        self.d_k = np.zeros_like(self.d_real, dtype=complex)
        self.d_real_is_computed = False
        self.d_k_is_computed = False

        # T-matrix in real and k-space
        self.t_real = np.zeros([3 * n_bas, 3 * n_bas, lat.l[0], lat.l[1], lat.l[2]])
        self.t_k = np.zeros_like(self.d_real, dtype=complex)
        self.t_real_is_computed = False
        self.t_k_is_computed = False

        # lattice sum energy
        self.egy = 0.0
        self.energy_is_computed = False

        # lattice sum pressure
        self.prs = 0.0
        self.prs_is_computed = False

        # running coefficient :math:`d ^{(1)}_s`
        self.coeff_d1 = np.zeros([self.lcoord[-1]])

        # running coefficient :math:`d^{(2)}_{s|\\alpha, \\beta}`
        self.coeff_d2 = np.zeros([3, 3, self.lcoord[-1]])

        mat_h = np.zeros([3, lat.l[0], lat.l[1], lat.l[2]])
        mat_b = np.zeros([3, n_bas, n_bas])

        # iterate over species
        self.b_grid = [(x, y) for x in range(n_spc) for y in range(n_spc)]
        # iterate over dimension
        self.a_grid = [(x, y) for x in range(3) for y in range(3)]

        # define the iteration over elements
        mat3ca, mat3cb, mat3cc = np.ix_(np.arange(lat.l[0]), np.arange(lat.l[1]), np.arange(lat.l[2]))

        for ind in range(3):
            mat_h[ind] = lat.a_vector[0, ind] * mat3ca + lat.a_vector[1, ind] * mat3cb + lat.a_vector[2, ind] * mat3cc

        # define primitive vectors
        for l1, l2 in self.b_grid:
            for b1 in range(lat.typ[l1]):
                for b2 in range(lat.typ[l2]):
                    ind1 = b1 + self.ncoord[l1]
                    ind2 = b2 + self.ncoord[l2]
                    self.r_vec[:, ind1, ind2, :, :, :] += mat_h[:, :, :, :]
        # define basis
        for l1, l2 in self.b_grid:
            for ind0 in range(3):
                for b1 in range(lat.typ[l1]):
                    for b2 in range(lat.typ[l2]):
                        ind1 = b1 + self.ncoord[l1]
                        ind2 = b2 + self.ncoord[l2]
                        mat_b[ind0, ind1, ind2] = lat.v_vector[ind1, ind0] - lat.v_vector[ind2, ind0]
        # compute positions by broadcasting the basis
        self.r_vec += mat_b[:, :, :, np.newaxis, np.newaxis, np.newaxis]

    def latsum_scalar(self, fun):
        """Returns the lattice sum for a given function

        :param fun: the given function
        :return: lattice sum
        :rtype: float
        """

        val = 0.0

        for l1, l2 in self.b_grid:
            b_mat = self.r_vec[:, self.ncoord[l1]:self.lcoord[l1], self.ncoord[l2]:self.lcoord[l2], :, :, :]
            val += np.sum(fun(self.d_d.min_dist(b_mat), (l1, l2)))

        return 0.5*val/self.num_species

    def coeff_linear(self):
        """Computes the linear coefficient (should be zero) in the expansion of the D-matrix

        :return: The linear coefficient
        :rtype: numpy.ndarray
        """

        alst = (1, 2, 3, 4)

        val = np.zeros([3, self.lcoord[-1]])

        for l1, l2 in self.b_grid:
            c_ind = np.arange(self.ncoord[l1], self.lcoord[l1])
            b_mat = self.r_vec[:, c_ind, self.ncoord[l2]:self.lcoord[l2], :, :, :]
            b_mat = self.d_d.min_dist(b_mat)
            for ind in range(3):
                val[ind, c_ind] += np.sum(self.pot.g1(b_mat, (l1, l2)) * b_mat[ind], axis=alst)

        return val

    def latsum_coef1(self, fun):
        """Lattice sum defining running coefficient :math:`d ^{(1)}_s` for an arbitrary function

        :param fun: the given function
        :return: the lattice sum
        :rtype: numpy.ndarray
        """

        val = np.zeros([self.lcoord[-1]])

        alst = (1, 2, 3, 4)

        for l1, l2 in self.b_grid:
            c_ind = np.arange(self.ncoord[l1], self.lcoord[l1])
            b_mat = self.d_d.min_dist(self.r_vec[:, c_ind, self.ncoord[l2]:self.lcoord[l2], :, :, :])
            val[c_ind] += np.sum(fun(b_mat, (l1, l2)), axis=alst)

        return val

    def latsum_coef2(self, fun):
        """Lattice sum defining running coefficient :math:`d^{(2)}_{s|\\alpha, \\beta}` for an arbitrary function

        :param fun: the given function
        :return: the lattice sum
        :rtype: numpy.ndarray
        """

        val = np.zeros([3, 3, self.lcoord[-1]])

        alst = (1, 2, 3, 4)

        for m1, m2 in self.a_grid:
            for l1, l2 in self.b_grid:
                c_ind = np.arange(self.ncoord[l1], self.lcoord[l1])
                b_mat = self.d_d.min_dist(self.r_vec[:, c_ind, self.ncoord[l2]:self.lcoord[l2], :, :, :])
                val[m1, m2, c_ind] += np.sum(fun(b_mat, (l1, l2)) * b_mat[m1] * b_mat[m2], axis=alst)

        return val

    def dmat_arbit(self, fun1, fun2):
        """Computes a D-matrix from two arbitrary functions

        :param fun1: first given function
        :param fun2: second given function
        :return: d-matrix (in real space)
        :rtype: numpy.ndarray
        """

        darb_c1 = self.latsum_coef1(fun1)
        darb_c2 = self.latsum_coef2(fun2)

        darb = np.zeros_like(self.d_real)

        for l1, l2 in self.b_grid:
                b_range = np.arange(self.ncoord[l1], self.lcoord[l1])
                c_ind1, c_ind2 = np.ix_(b_range, np.arange(self.ncoord[l2], self.lcoord[l2]))
                b_mat = self.d_d.min_dist(self.r_vec[:, c_ind1, c_ind2, :, :, :])
                for m1 in range(3):
                    darb[m1 + 3*c_ind1, m1 + 3*c_ind2] -= fun1(b_mat, (l1, l2))
                    for m2 in range(3):
                        darb[m1 + 3*c_ind1, m2 + 3*c_ind2] -= fun2(b_mat, (l1, l2))*b_mat[m1]*b_mat[m2]

        # add the zero mode
        for l1 in range(self.n_spc):
            c_ind = np.arange(self.ncoord[l1], self.lcoord[l1])
            for m1 in range(3):
                darb[m1 + 3*c_ind, m1 + 3*c_ind, 0, 0, 0] += darb_c1[c_ind]
                for m2 in range(3):
                    darb[m1 + 3*c_ind, m2 + 3*c_ind, 0, 0, 0] += darb_c2[m1, m2, c_ind]

        return [darb_c1, darb_c2, darb]

    def energy(self):
        """Returns the value of the energy (lattice sum)

        :return: value of the lattice sum energy
        :rtype: float
        """

        if self.energy_is_computed:
            return self.egy

        self.egy = self.latsum_scalar(self.pot.e1)

        self.energy_is_computed = True

        return self.egy

    def dmatrix_real(self):
        """Computes the D-matrix in real space

        :return: The D-matrix in real space
        :rtype: numpy.ndarray
        """

        # avoid computing the D-matrix twice
        if self.d_real_is_computed:
            return self.d_real

        val = self.dmat_arbit(self.pot.g1, self.pot.g2_over_r2)

        self.coeff_d1 = val[0]
        self.coeff_d2 = val[1]

        self.d_real = val[2]

        self.d_real_is_computed = True

        return self.d_real

    def prs_enthalpic(self):
        """Returns the value of :math:`\\frac{P_0 V}{N}` (lattice sum)

        :return: value of the lattice sum energy
        :rtype: float
        """

        if self.prs_is_computed:
            return self.prs

        self.prs = self.latsum_scalar(self.pot.e2)/3.0

        self.prs_is_computed = True

        return self.prs

    def dmatrix_k(self):
        """Computes the D-matrix in k-space

        :return: The D-matrix in k-space
        :rtype: numpy.ndarray
        """

        if self.d_k_is_computed:
            return self.d_k

        self.d_k = fftn(self.dmatrix_real(), axes=(2, 3, 4))

        self.d_k_is_computed = True

        return self.d_k

    def tmatrix_real(self):
        """Computes the T-matrix in real space

        :return: T-matrix in real space
        :rtype: numpy.ndarray
        """
        # avoid computing the T-matrix twice
        if self.t_real_is_computed:
            return self.t_real

        val = self.dmat_arbit(self.pot.g2, self.pot.g3_over_r2)

        self.t_real = val[2]

        self.t_real_is_computed = True

        return self.t_real

    def tmatrix_k(self):
        """Computes the T-matrix in k-space

        :return: The T-matrix in k-space
        :rtype: numpy.ndarray
        """

        if self.t_k_is_computed:
            return self.t_k

        self.t_k = fftn(self.tmatrix_real(), axes=(2, 3, 4))

        self.t_k_is_computed = True

        return self.t_k
