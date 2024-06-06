"""
:module: Quaternion
:platform: Unix, Windows
:synopsis: Quaternion class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  May 2017
"""

from __future__ import division
import numpy as np


class Quaternion(object):
    """
    Defines the class for quaternion including all of their operations
    """

    def __init__(self, val):
        """The constructor
        
        :param val: Mx4 dimensional numpy array
        """

        # make sure that the matrix is of type (M, 4)
        if len(val.shape) == 1:
            self.quat = np.zeros([1, 4])
            self.quat[0, :] = val
        else:
            self.quat = val
        self.dim = val.shape[0:-1]

    def __add__(self, other):
        """Quaternion addition
        
        :param other: Quaternion object
        :return: addition
        :rtype: Quaternion
        """

        qsum = Quaternion(np.zeros_like(self.quat))
        qsum.quat = self.quat + other.quat

        return qsum

    def __str__(self):
        """String representation
        """
        return np.array_str(self.quat)

    def __neg__(self):
        """
        Negative value
        :return: - value
        :rtype: Quaternion
        """
        qneg = Quaternion(np.zeros_like(self.quat))
        qneg.quat = - self.quat

        return qneg

    def __eq__(self, other):
        """Defines equality
        
        :param other: Quaternion object
        :return: True if both quaternions are the same
        :rtype: bool
        """
        return np.all(np.equal(self.quat, other.quat))

    def __mul__(self, other):
        """Quaternion right multiplication
        
        :param other: Quaternion object
        :return: the quantity q*other 
        :rtype: quaternion
        """

        mult = Quaternion(np.zeros_like(self.quat))

        mult.quat[..., 0] = self.quat[..., 0] * other.quat[..., 0] - self.quat[..., 1] * other.quat[..., 1] - \
                            self.quat[..., 2] * other.quat[..., 2] - self.quat[..., 3] * other.quat[..., 3]

        mult.quat[..., 1] = self.quat[..., 0] * other.quat[..., 1] + self.quat[..., 1] * other.quat[..., 0] + \
                            self.quat[..., 2] * other.quat[..., 3] - self.quat[..., 3] * other.quat[..., 2]

        mult.quat[..., 2] = self.quat[..., 0] * other.quat[..., 2] - self.quat[..., 1] * other.quat[..., 3] + \
                            self.quat[..., 2] * other.quat[..., 0] + self.quat[..., 3] * other.quat[..., 1]

        mult.quat[..., 3] = self.quat[..., 0] * other.quat[..., 3] + self.quat[..., 1] * other.quat[..., 2] - \
                            self.quat[..., 2] * other.quat[..., 1] + self.quat[..., 3] * other.quat[..., 0]

        return mult

    def __sub__(self, other):
        """Quaternion subtraction
        
        :param other: Quaternion object
        :return: Sutraction
        :rtype: Quaternion
        """
        qdiff = self + (-other)

        return qdiff

    def __truediv__(self, other):
        """Quaternion division

        :param other: Quaternion object
        :return: Division
        :rtype: Quaternion
        """

        divn = self * other.inverse()

        return divn

    def __div__(self, other):
        """Quaternion division
        
        :param other: Quaternion object
        :return: Division
        :rtype: Quaternion
        """

        divn = self*other.inverse()

        return divn

    def conjugate(self):
        """Conjugate value
        
        :return: The conjugation
        :rtype: Quaternion
        """

        conj = Quaternion(np.zeros_like(self.quat))
        conj.quat[..., 0] = self.quat[..., 0]
        conj.quat[..., 1:] = - self.quat[..., 1:]

        return conj

    def norm(self):
        """Norm of the quaternion
        
        :return: The quaternion norm
        :rtype: ndarray
        """

        return np.sqrt(np.sum(self.quat**2, axis=-1))

    def inverse(self):
        """Inverse of a quaternion
        
        :return: The inverse
        :rtype: Quaternion
        """
        inv = Quaternion(np.zeros_like(self.quat))
        norm = np.zeros_like(self.quat)

        for ind in range(4):
            norm[..., ind] = self.norm()

        inv.quat = self.conjugate().quat/norm**2
        return inv

    def approx(self, other, eps=1e-10):
        """Defines whether two quaternions are close
        
        :param other: Quaternion object
        :param eps: closeness approach
        :return: Whether the condition is true
        :rtype: bool
        """

        return np.allclose(self.quat, other.quat, eps)


def rotation_matrix(quaternion):
    """Gives rotation matrix for given quaternion
    :param quaternion:
    :return: rotation matrix
    """
    a, b, c, d = quaternion
    mat = np.array([[a**2+b**2-c**2-d**2, 2*(b*c+a*d), 2*(b*d-a*c)],
                    [2*(b*c-a*d), a**2-b**2+c**2-d**2, 2*(c*d+a*b)],
                    [2*(b*d+a*c), 2*(c*d-a*b), a**2-b**2-c**2+d**2]]).T
    return mat
