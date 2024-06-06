"""
:module: ThermFunctions
:platform: Unix, Windows
:synopsis: Computes thermodynamic functions that do not have explicit temperature dependence

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, January2016
"""

import numpy as np

import hoodlt.D_matrix.Dmatrix_Mixt_D as Dm
import numpy.linalg as la
from numpy.fft import ifftn


class ThermoFunction(object):
    """
        Defines the class that computes thermodynamic functions that do not have explicit temperature dependence
    """

    def __init__(self, lat, pot, res=1e-8):
        """The constructor

        :param lat: lattice object
        :param pot: potential object
        :param res: eigenvalue resolution
        """

        self.pot = pot
        self.lat = lat
        self.num_pnts = float(self.lat.num_pnts())
        self.d_mat = Dm.DMatrix(lat, pot)

        self.egy = self.d_mat.energy()
        self.prs_ent = self.d_mat.prs_enthalpic()
        self.epy = None
        self.prs_epy = None

        self.mat_d = self.d_mat.dmatrix_k()
        self.mat_t = self.d_mat.tmatrix_k()

        self.inv_d = np.zeros_like(self.mat_d, dtype=complex)

        mat_indx = self.mat_d.shape
        num_eigvals = mat_indx[1]*mat_indx[2]*mat_indx[3]*mat_indx[4]
        self.eigvals = np.zeros([mat_indx[2], mat_indx[3], mat_indx[4], mat_indx[1]])

        for ind1 in range(mat_indx[2]):
            for ind2 in range(mat_indx[3]):
                for ind3 in range(mat_indx[4]):
                    self.eigvals[ind1, ind2, ind3] = la.eigvalsh(self.mat_d[:, :, ind1, ind2, ind3])

        # count the number of positive eigenvalues
        num_positive = np.sum(self.eigvals > res)

        # count the number of negative eigenvalues
        num_negative = np.sum(self.eigvals < -res)

        # count the number of zero eigenvalues
        num_zero = num_eigvals - num_positive - num_negative

        self.eigenvalues_type = np.array([num_positive, num_zero, num_negative], dtype=int)

        self.stable = False
        if (self.eigenvalues_type[1] == 3) and (self.eigenvalues_type[2] == 0):
            self.stable = True

        if self.stable:
            self.epy = -0.5*np.sum(np.log(self.eigvals[self.eigvals > res]))/self.num_pnts

            for ind1 in range(mat_indx[2]):
                for ind2 in range(mat_indx[3]):
                    for ind3 in range(mat_indx[4]):
                        if not (ind1 == 0 and ind2 == 0 and ind3 == 0):
                            self.inv_d[:, :, ind1, ind2, ind3] = la.inv(self.mat_d[:, :, ind1, ind2, ind3])

            veig, weig = la.eigh(self.mat_d[:, :, 0, 0, 0])
            non_zero = np.where(veig > res)

            for ind_v in non_zero[0]:
                self.inv_d[:, :, 0, 0, 0] += np.outer(np.conjugate(weig)[:, ind_v], weig[:, ind_v]/veig[ind_v])

            val = 0.0
            for ind1 in range(mat_indx[2]):
                for ind2 in range(mat_indx[3]):
                    for ind3 in range(mat_indx[4]):
                        amat = self.inv_d[:, :, ind1, ind2, ind3]
                        bmat = self.mat_t[:, :, ind1, ind2, ind3]
                        val += np.trace(np.dot(amat, bmat))

            self.prs_epy = -np.real(val)/(6*self.num_pnts)

            self.inv_d_real = np.real(ifftn(self.inv_d, axes=(2, 3, 4)))

    def energy(self):
        """Returns the lattice sum energy

        :return: energy per particle (in the same units as the potential is defined)
        :rtype: float
        """

        return self.egy

    def prs_enthalpic(self):
        """Returns the quantity :math:`\\frac{P_0 V}{N}`

        :return: the pressure (in energy units, units given by the potential)
        :rtype: float
        """

        return self.prs_ent

    def entropy(self):
        """Returns the entropy

        :return: Entropy in units of :math:`k_B`
        :rtype: float
        """

        return self.epy

    def prs_entropic(self):
        """Returns the quantity :math:`\\frac{P_D V}{N}`

        :return: the pressure (in energy units, units given by the potential)
        :rtype: float
        """

        return self.prs_epy

    def prs_total(self):
        """Returns the total pressure

        :return: :math:`\\frac{(P_0+P_D)}{N V}`
        :rtype: float
        """

        return self.prs_ent + self.prs_epy

    def dmat_vib(self, nmat):
        """Vibration matrix: Computes the average displacement

        :param nmat: integers describing the elements
        :return: vibration matrix
        :rtype: numpy.ndarray
        """

        return self.inv_d_real[:, :, nmat[0], nmat[1], nmat[2]]
