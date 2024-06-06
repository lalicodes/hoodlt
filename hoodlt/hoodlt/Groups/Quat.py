"""
:module: QuaternionRowan
:platform: Unix, Windows
:synopsis: Quaternion class

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  June 2022
"""


import rowan as rw


class Quat(object):
    """
    Defines the class for quaternion including all of their operations
    """

    def __init__(self):
        """The constructor """

    @staticmethod
    def rotation(qt, vec):
        """
        :param qt: quaternion as a tuple
        :param vec: numpy array

        :return : rotated vectors as numpy array
        """

        return rw.rotate(qt, vec)

    @staticmethod
    def vector_vector_rotation(v1, v2):
        """
        :param v1: vector the quaternion goes from
        :param v2: vector the quaternion goes to

        :return : rotated vectors as numpy array
        """

        return rw.vector_vector_rotation(v1, v2)

    @staticmethod
    def multiply(v1, v2):
        """
        :param v1: vector the quaternion goes from
        :param v2: vector the quaternion goes to

        :return : quaternion multiplication
        """

        return rw.multiply(v1, v2)

    @staticmethod
    def inverse(v1):
        """
        :param v1: vector the quaternion goes from

        :return : quaternion inverse
        """

        return rw.inverse(v1)

    @staticmethod
    def from_axis_angle(axes, angles):
        """
        :param axes: An array of vectors (the axes of rotation)
        :param angles: An array of angles in radians, assumed counterclockwise about the axes

        :return : rotations as quaternions in array
        """

        return rw.from_axis_angle(axes, angles)

    @staticmethod
    def single_from_axis_angle(axes, angles):
        """
        :param axes: An array of vectors (the axes of rotation)
        :param angles: An array of angles in radians, assumed counterclockwise about the axes

        :return : rotations as a single quaternion, applying the rotations in the order given
        """
        
        all_quats = list(rw.from_axis_angle(axes, angles))
        total_quat = all_quats.pop(0)
        for quat in all_quats:
            total_quat = rw.multiply(quat, total_quat)

        return total_quat
