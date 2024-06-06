"""
:module: Dmatrix_Mixt_potential_simple
:platform: Unix, Windows
:synopsis: Defines the abstract class for the lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
"""

import abc
from hoodlt.D_matrix.Dmatrix_Mixt_Potential_Decorator import compute_norm_keep_matrix


class DPotential(object):
    """Defines a general pair potential
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        """The constructor"""
        self.cut = 0.0
        self.epsilon = 0.0
        self.indx = (0, 0)

    def __repr__(self):
        """Potential name

        :return: potential name
        :rtype: str
        """
        return 'Potential name'
        
    @abc.abstractmethod
    def der0(self, r):
        """Potential value

        :param r: point
        :return: value
        :rtype: float
        """

        return r
  
    @abc.abstractmethod
    def der1(self, r):
        """Potential 1st-derivative value

        :param r: point
        :return: value
        :rtype: float
        """
        return r

    @abc.abstractmethod
    def der2(self, r):
        """Potential 2nd-derivative value

        :param r: point
        :return: value
        :rtype: numpy.ndarray
        """
        return r

    @abc.abstractmethod
    def der3(self, r):
        """Potential 3rd-derivative value

        :param r: point
        :return: value
        :rtype: float
        """
        return r

    @compute_norm_keep_matrix
    def e1(self, *args):
        """The function :math:`e_1(r)` necessary to compute energy

        :param args: arguments
        :return: energy function
        """
        return self.der0(args[0])

    @compute_norm_keep_matrix
    def e2(self, *args):
        """The function :math:`e_2(r)` necessary to compute pressure

        :param args: arguments
        :return: enthalpy pressure function
        """
        return -args[0]*self.der1(args[0])

    @compute_norm_keep_matrix
    def g1(self, *args):
        """The :math:``g_1(r)`` value

        :param args: arguments
        :return: energy function
        """

        return self.der1(args[0])/args[0]

    @compute_norm_keep_matrix
    def g2(self, *args):
        """The :math:``g_2(r)`` value

        :param args: arguments
        :return: energy function
        """

        return self.der2(args[0]) - self.der1(args[0])/args[0]

    @compute_norm_keep_matrix
    def g2_over_r2(self, *args):
        """The :math:``g_2(r)/r^2`` value

        :param args: arguments
        :return: energy function
        """
        r = args[0]
        return self.der2(r)/r**2 - self.der1(r)/r**3

    @compute_norm_keep_matrix
    def g3(self, *args):
        """The :math:``g_3(r)`` value

        :param args: arguments
        :return: energy function
        """
        r = args[0]
        return r*self.der3(r) - self.der2(r) + self.der1(r)/r

    @compute_norm_keep_matrix
    def g3_over_r2(self, *args):
        """The :math:``g_3(r)//r^2`` value

        :param args: arguments
        :return: energy function
        """
        r = args[0]
        return self.der3(r)/r - self.der2(r)/r**2 + self.der1(r)/r**3
