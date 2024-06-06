"""
:module: D_matrix_Mixt_Lattices
:platform: Unix, Windows
:synopsis: Defines the abstract class for the lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""

import numpy as np
import abc


class DMixtLattice(object):

    """Defines a lattice (this class is abstract)

    The lattice is defined by
    
    .. math::
    
      {\\bm R}_{(n,s)}={\\bm v}_{s}+\sum_{l=1}^3 n_l {\\bm a}_l

    where
    
    .. math::
    
      {\\bm v}_{s} \\quad s=1\\cdots m

    defines the basis
    
    .. math::
    
      n_l = 0 \\cdots L-1

    defines the number of unit cells
    
    .. math::
    
     {\\bm a}_l \\quad l= 1,2,3

    defines the primitive vectors
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, l_value, a_nn_ext):
        """The Constructor

        :param l_value: linear size
        :param a_nn_ext: effective diameter of A particle
        :return:
        """

        if isinstance(l_value, list):
            [l1, l2, l3] = l_value
        elif isinstance(l_value, int):
            # Box length along the 1st primitive vector
            l1 = l_value
            # Box length along the 2nd primitive vector
            l2 = l_value
            # Box length along the 3rd primitive vector
            l3 = l_value
        else:
            raise ValueError('first argument of constructor is not correct')

        self.l = np.array([l1, l2, l3])

        self.a_nn = a_nn_ext

        self.a_axis = self.a_nn
        self.b_axis = self.a_nn
        self.c_axis = self.a_nn

        self.a_vector = np.zeros((3, 3))
        self.v_vector = np.zeros((1, 3))
        self.typ = np.zeros((1,), dtype=int)
        self.radius = np.ones((1,))

    def a_val(self):
        """Minimum distance of hard sphere A-A particles

        :return: minimal distance
        :rtype: float
        """
        return self.a_nn

    def get_a(self, l):
        """Primitive vector of corresponding index

        :param l: base index
        :return: base corresponding to l
        :rtype: numpy.ndarray
        """
        return self.a_vector[l]

    def get_v(self, k):
        """Basis vector (matrix)

        :param k: index of the corresponding base
        :return: corresponding base
        :rtype: numpy.ndarray
        """

        return self.v_vector[k[1]+np.sum(self.typ[0:k[0]])]

    def num_pnts(self):
        """ Number of particles within the box

        :return: number of particles as [int num A,int num B, int num C, etc..]
        :rtype: list
        """
        return np.prod(self.l)*self.typ

    def l_box(self):
        """ Returns the box size in the natural unit of length

        :return: lattice box in units
        :rtype: numpy.ndarray
        """
        return np.array([self.a_axis, self.b_axis, self.c_axis])*self.l

    def vol_unit_cell(self):
        """Volume of the unit cell thus defined

        :return: Volume of the unit cell
        :rtype: float
        """

        return np.dot(np.cross(self.get_a(0), self.get_a(1)), self.get_a(2))

    def redefine_lattice_constant(self, a_lat):
        """Allows to redefine the lattice constant for the lattice, without changing the a_nn parameter, that is, the
        units of the lattice constant.

        :param a_lat: new lattice constant, it can be a scalar or a numpy array, which allows to scale the three axis
        """

        if isinstance(a_lat, np.ndarray):
            vec = a_lat
        elif isinstance(a_lat, list):
            vec = np.array(a_lat)
        else:
            vec = np.array([a_lat, a_lat, a_lat])

        fac = vec/(self.a_nn*np.array([self.a_axis, self.b_axis, self.c_axis]))

        self.a_axis *= fac[0]
        self.b_axis *= fac[1]
        self.c_axis *= fac[2]

        self.a_vector = np.transpose(np.transpose(self.a_vector)*fac)
        self.v_vector = self.v_vector*fac

    @staticmethod
    def name():
        """Lattice full name

        :return: [name, latex name, space group, strukturberich]
        :rtype: list
        """
        return ['name', 'latex_name', 'space_group', 'strukturberich']

    def g_l(self):
        """ Density function, the density :math:`\\frac{N}{V} = \\frac{g_l}{a^3_{nn}}`

        :return: density coefficient
        :rtype: float
        """
        return np.sum(self.typ)*self.a_nn**3/self.vol_unit_cell()

    @abc.abstractmethod
    def pf(self):
        """Packing fraction

        :return: packing fraction of the corresponding lattice
        :rtype: float
        """
        return None
