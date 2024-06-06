"""
:module: LatticePoints
:platform: Unix, Windows
:synopsis: Class that labels points in a lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> February 2021
"""

import numpy as np
import numpy.linalg as la


class LatPoints(object):
    """ Defines the class identifying each lattice point by either one integer :math:`i_{val}` or by the more
    intuitive five integers :math:`(i_b,j_b,l_1,l_2,l_3)` as described in documentation
    """

    def __init__(self, lat):
        """The constructor

               :param lat: lattice object
        """
        self.lat = lat

        # number of lattice species
        self.n_spc = len(lat.typ)
        # number of basis in the unit cell
        self.n_bas = np.sum(lat.typ)

        # define internal coordinates
        self.ncoord = np.zeros(self.n_spc, dtype=int)
        self.lcoord = np.zeros(self.n_spc, dtype=int)
        bfac = 0
        for dims in range(self.n_spc):
            self.ncoord[dims] = bfac
            bfac += lat.typ[dims]
            self.lcoord[dims] = bfac

    def i_val(self, ind_pnt):
        """Returns the single coordinate site for the lattice

            :param ind_pnt: ndarray containing 5 integers :math:`(i_b,j_b,l_1,l_2,l_3)`
            :return: single variable of a point (:math:`i_{val}`)
        """
        vec = np.array([self.ncoord[ind_pnt[0]], 1, 1, self.lat.l[0], self.lat.l[0]*self.lat.l[1]])
        vec_b = np.array([1, 1, self.lcoord[-1], self.lcoord[-1], self.lcoord[-1]])
        return np.sum(vec*vec_b*ind_pnt)

    def coordinate_five(self, i_pnt):
        """Returns the five coordinate label for the lattice

            :param i_pnt: integer, variable :math:`i_{val}`
            :return: n_darray 5 integers :math:`(i_b,j_b,l_1,l_2,l_3)`
        """

        coor = np.zeros(5, dtype='int')

        coor[-1] = i_pnt//(self.n_bas*self.lat.l[0]*self.lat.l[1])
        coor[-2] = (i_pnt-coor[-1]*self.n_bas*self.lat.l[0]*self.lat.l[1]) // (self.n_bas * self.lat.l[0])
        coor[-3] = (i_pnt-self.n_bas*self.lat.l[0]*(coor[-2]+self.lat.l[1] * coor[-1]))//self.n_bas
        num = i_pnt - self.n_bas*(coor[-3]+coor[-2]*self.lat.l[0]+coor[-1]*self.lat.l[0]*self.lat.l[1])

        index_i = np.argwhere(np.logical_and((num < self.lcoord), num >= self.ncoord))

        coor[0] = index_i[0]
        coor[1] = num - self.ncoord[coor[0]]

        return coor
