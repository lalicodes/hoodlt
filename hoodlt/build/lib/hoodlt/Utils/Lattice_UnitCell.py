# noinspection PyInterpreter
"""
:module: Lattice_UnitCell
:platform: Unix, Windows
:synopsis: Utility to isolate the geometric properties of a unit cell

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2015
"""

import numpy as np
import hoodlt.Utils.LatticeCoords as Lc


class UnitCell(object):
    """Unit cell class: it defines
    
    :ivar pnt: numpy array containing points within the unit cell and nearest neighbors
    :ivar indices_all: type of vector
    :ivar u_cell: indices of the unit cell
    :ivar basis_coord: position of point given the coordinate basis
    :ivar lat: lattice object
    :ivar ratio_b: relative radius of B particles :math:`\\frac{\\bar{r}_B}{r_B}`
    """

    def __init__(self, lat):
        """The constructor
        
        :param lat: :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
        """

        # lattice
        self.lat = lat
        ntotal = 27 * np.sum(lat.typ)
        # total number of points, the ones defining the unit cell and additional points including all nearest neighbors
        self.num_points = ntotal
        # Store all the points
        self.pnts = np.zeros([ntotal, 3])
        # indices of all points
        self.indices_all = np.zeros([ntotal, 2], dtype=int)
        # number of basis (number of points within the unit cell)
        self.num_basis = np.sum(lat.typ)
        # indices of the basis
        self.u_cell = np.zeros([self.num_basis, 3], dtype=int)
        # position of unit cell points from index (p, q)
        self.basis_coord = [np.zeros(vl, dtype=int) for vl in lat.typ]

        self.frac = Lc.LatCords(self.lat).fractional

        n_grid = [(n0, n1, n2) for n0 in range(-1, 2, 1) for n1 in range(-1, 2, 1) for n2 in range(-1, 2, 1)]

        ind_total = 0
        ind_ucell = 0

        for ind_1 in range(len(lat.typ)):
            for ind_2 in range(lat.typ[ind_1]):
                for n0, n1, n2 in n_grid:
                    pn = lat.get_v([ind_1, ind_2]) + n0 * lat.get_a(0) + n1 * lat.get_a(1) + n2 * lat.get_a(2)
                    self.pnts[ind_total, :] = pn[:]
                    self.indices_all[ind_total, 0] = ind_1
                    self.indices_all[ind_total, 1] = ind_2
                    xv = self.frac(pn)
                    if ((xv[0] < 1) and xv[0] >= 0) and ((xv[1] < 1) and xv[1] >= 0) and ((xv[2] < 1) and xv[2] >= 0):
                        self.u_cell[ind_ucell, 0] = ind_total
                        self.u_cell[ind_ucell, 1:] = self.indices_all[ind_total, :]
                        self.basis_coord[ind_1][ind_2] = ind_total
                        ind_ucell += 1
                    ind_total += 1

    def bounding_box(self):
        """
        Computes a bounding box such that all points are contained within
        
        :return:  max and minimum along each direction
        :rtype: numpy array
        """

        # set up a bounding box
        val1 = 3 * self.lat.get_a(0) + 3 * self.lat.get_a(1) + 3 * self.lat.get_a(2)

        mat = np.array([val1, -val1])

        return np.array([np.amin(mat, axis=0), np.amax(mat, axis=0)])


def unit_cell(lat, max_ind=1):
    """Returns a unit cell structure

    The unit cell structure is a dictionary, with the following entries:
    
    points: all points within the unit cell and their neighbors
    
    indices_all: particle type, index within the type
    
    indices: it is a matrix with three entries: Particle index of the unit cell, particle type and number within type
    
    lattice: lattice class

    :param lat: :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
    :param max_ind: maximum_indices to return
    :return: points within unit cell and all its neighbors, [Particle index, particle type, particle number, lattice]
    :rtype: dict
    """

    n_grid = [(n0, n1, n2) for n0 in range(-max_ind, max_ind+1, 1) for n1 in range(-max_ind, max_ind+1, 1)
              for n2 in range(-max_ind, max_ind+1, 1)]

    pnts = np.zeros([((2*max_ind+1)**3)*np.sum(lat.typ), 3])
    index_all = []

    pnt_unit_cell = []

    ind_total = 0
    for ind_1 in range(len(lat.typ)):
        for ind_2 in range(lat.typ[ind_1]):
            for n0, n1, n2 in n_grid:
                r_v = lat.get_v([ind_1, ind_2]) + n0*lat.get_a(0) + n1*lat.get_a(1) + n2*lat.get_a(2)
                pnts[ind_total, :] = r_v[:]
                index_all.append([ind_1, ind_2, n0, n1, n2])
                
                if (n0 == 0) and (n1 == 0) and (n2 == 0):
                    pnt_unit_cell.append([ind_total, ind_1, ind_2])
                    
                ind_total += 1

    return {'points': pnts, 'indices_all': index_all, 'indices': pnt_unit_cell, 'lattice': lat}
