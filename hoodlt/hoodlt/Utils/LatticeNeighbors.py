"""
:module: LatticeNeighbors
:platform: Unix, Windows
:synopsis: Class that defines positions, types and nearest neighbors for any lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> January 2020
"""

import numpy as np
import numpy.linalg as la
import hoodlt.D_matrix.Dmatrix_Mixt_distances as Dm


class LatNeighbor(object):
    """ Defines the class providing lattice sites and  nearest neighbors. It may be used to build lattices in hoomd
    and also add springs to enhance the lattice stability
    """

    def __init__(self, lat):
        """The constructor

               :param lat: lattice object
        """

        self.lat = lat
        # cut-off is slightly less (7/10) than L/2
        cut_off = np.amin(7*lat.a_nn*lat.l)/20
        # define distances with the minimum convention
        self.d_d = Dm.MinConv(lat, cut_off)
        # number of digits (necessary for np.round)
        self.epsl = 8

        # number of lattice species
        n_spc = len(lat.typ)
        # number of basis in the unit cell
        n_bas = np.sum(lat.typ)
        # store distances for each base in a list
        self.u_dist = n_bas*[0]

        # define internal coordinates
        self.ncoord = np.zeros(n_spc, dtype=int)
        self.lcoord = np.zeros(n_spc, dtype=int)
        bfac = 0
        for dims in range(n_spc):
            self.ncoord[dims] = bfac
            bfac += lat.typ[dims]
            self.lcoord[dims] = bfac

        # point positions
        self.l_vec = np.zeros([3, np.sum(lat.num_pnts())])
        r_vec = np.zeros([3, n_bas, n_bas, lat.l[0], lat.l[1], lat.l[2]])
        # point types and base types
        self.typ = np.zeros([np.sum(lat.num_pnts())], dtype='int')
        self.base = np.zeros([np.sum(lat.num_pnts())], dtype='int')

        # helper matrices
        mat_h = np.zeros([3, lat.l[0], lat.l[1], lat.l[2]])
        mat_b = np.zeros([3, n_bas, n_bas])

        # define the iteration over elements
        mat3ca, mat3cb, mat3cc = np.ix_(np.arange(lat.l[0]), np.arange(lat.l[1]), np.arange(lat.l[2]))

        # define the translations of the unit cell
        for ind in range(3):
            mat_h[ind] = lat.a_vector[0, ind] * mat3ca + lat.a_vector[1, ind] * mat3cb + lat.a_vector[2, ind] * mat3cc

        # iterate over species
        self.b_grid = [(x, y) for x in range(n_spc) for y in range(n_spc)]
        # iterate over dimension
        self.a_grid = [(x, y) for x in range(3) for y in range(3)]

        # define primitive vectors
        for l1, l2 in self.b_grid:
            for b1 in range(lat.typ[l1]):
                for b2 in range(lat.typ[l2]):
                    ind1 = b1 + self.ncoord[l1]
                    ind2 = b2 + self.ncoord[l2]
                    r_vec[:, ind1, ind2, :, :, :] += mat_h[:, :, :, :]
        # define basis
        for l1, l2 in self.b_grid:
            for ind0 in range(3):
                for b1 in range(lat.typ[l1]):
                    for b2 in range(lat.typ[l2]):
                        ind1 = b1 + self.ncoord[l1]
                        ind2 = b2 + self.ncoord[l2]
                        mat_b[ind0, ind1, ind2] = lat.v_vector[ind1, ind0] - lat.v_vector[ind2, ind0]
        # compute positions by broadcasting the basis
        r_vec += mat_b[:, :, :, np.newaxis, np.newaxis, np.newaxis]

        self.dist = la.norm(self.d_d.min_dist(r_vec), axis=0)

        # sorted list with all distances for each representative basis point
        for ind in range(n_bas):
            self.u_dist[ind] = np.sort(np.unique(np.round(self.dist[ind], decimals=self.epsl)))

        mat_h = np.reshape(mat_h, (3, np.prod(lat.l)))
        mat_ind = np.arange(np.prod(lat.l))
        # create points
        ind_b = 0
        for ind_s in range(n_spc):
            for ind_t in range(lat.typ[ind_s]):
                new_ind = ind_b + n_bas*mat_ind
                vc = lat.get_v([ind_s, ind_t])
                self.l_vec[:, new_ind] = self.d_d.min_dist(vc[:, np.newaxis] + mat_h)
                self.typ[new_ind] = ind_s
                self.base[new_ind] = ind_b
                ind_b += 1

    def neighbor(self, ind_pnt, degree):
        """Returns the neighbor of degree(1=nearest neighbor,2=next to nearest neighbor, etc..) for site ind_pnt

            :param ind_pnt: index of a lattice point
            :param degree: degree is nearest neighbor=1, next to nearest neighbor = 2, etc..
            :return: ndarray with the neighbor indices (i_val)
        """

        if degree < 1:
            return []

        num_base = self.base[ind_pnt]
        dist = self.u_dist[num_base][degree]
        mat_c = dist*np.ones([self.l_vec.shape[1]])

        nv = self.l_vec[:, ind_pnt]

        nn = np.where(np.isclose(la.norm(self.d_d.min_dist(self.l_vec-nv[:, np.newaxis]), axis=0), mat_c))

        return nn[0]

    def bonds_in_lattice(self, p_type, degree):
        """
        Returns all pairs where base = ind_base forms bonds with degree neighbors

        :param p_type: base particle type (=0 for single component, 0,1 for binary, etc..)
        :param degree: degree  is nearest neighbor=1, next to nearest neighbor = 2, etc..
        :return: list of bonds, for example [[0,1], [2, 5], ..]
        :rtype: list
        """

        # list of all points with typ = p_type
        p_list = np.where(self.typ == p_type)[0]

        # list containing the bonds
        bond_list = []

        # build the list
        for ind_b in p_list:
            nn = self.neighbor(ind_b, degree)
            for ind_n in nn:
                if (self.typ[ind_n] == self.typ[ind_b]) and (ind_n < ind_b):
                    bond_list.append([ind_n, ind_b])
                elif self.typ[ind_n] != self.typ[ind_b]:
                    bond_list.append([ind_n, ind_b])

        return bond_list
