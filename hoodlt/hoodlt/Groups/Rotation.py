"""
:module: Rotation
:platform: Unix, Windows, OS
:synopsis: Defines the Polytope 335

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  February 2020
"""

from __future__ import division
import numpy as np
import hoodlt.Groups.Quaternion as Quat


class Rotate(object):
    """
     Defines rotations from Euler angles
    """

    def __init__(self, eul):
        """The Constructor
        
        :param eul: numpy array with the three Euler angles :math:`(\\phi, \\theta, \\psi)` as defined in Goldstein
        """

        mat_d = np.array([[np.cos(eul[0]), np.sin(eul[0]), 0], [-np.sin(eul[0]), np.cos(eul[0]), 0], [0, 0, 1]])
        mat_c = np.array([[1, 0, 0], [0, np.cos(eul[1]), np.sin(eul[1])], [0, -np.sin(eul[1]), np.cos(eul[1])]])
        mat_b = np.array([[np.cos(eul[2]), np.sin(eul[2]), 0], [-np.sin(eul[2]), np.cos(eul[2]), 0], [0, 0, 1]])

        self.space_to_body = np.dot(np.dot(mat_b, mat_c), mat_d)
        self.body_to_space = np.transpose(self.space_to_body)

        q0 = np.cos(0.5*eul[1])*np.cos(0.5*(eul[0]+eul[2]))
        q1 = np.sin(0.5*eul[1])*np.cos(0.5*(eul[0]-eul[2]))
        q2 = np.sin(0.5*eul[1])*np.sin(0.5*(eul[0]-eul[2]))
        q3 = np.cos(0.5*eul[1])*np.sin(0.5*(eul[0]+eul[2]))

        # quaternionic form
        self.quaternion = Quat.Quaternion(np.array([q0, q1, q2, q3]))

        # compute axis and angle of rotation
        cos_ang_half = q0
        sin_ang_half = np.sqrt(1-cos_ang_half**2)

        if np.abs(sin_ang_half) > 1e-8:
            self.axis_rotation = np.array([q1, q2, q3])/sin_ang_half
            self.angle_rotation = 2*np.arccos(cos_ang_half)
        else:
            self.axis_rotation = np.array([0.0, 0.0, 1.0])
            self.angle_rotation = 0.0


class RotateFromAxis(object):
    """
        Defines rotations from axis and angle
       """
    def __init__(self, naxis, ang):
        """
        Initialize the class from quaternions instead of euler angles

        :param naxis: ndarray axis of rotation
        :param ang: angle
        """

        self.axis_rotation = naxis
        self.angle_rotation = ang

        pr = np.sin(0.5*ang)
        q0 = np.cos(0.5*ang)
        q1 = pr*naxis[0]
        q2 = pr*naxis[1]
        q3 = pr*naxis[2]

        self.quaternion = Quat.Quaternion(np.array([q0, q1, q2, q3]))

        a1 = np.outer(naxis, naxis)
        a2 = np.eye(3)
        a3 = np.array([[0, naxis[2], -naxis[1]], [-naxis[2], 0, naxis[0]], [naxis[1], -naxis[0], 0]])

        self.space_to_body = a1 + (a2-a1)*np.cos(ang) + a3*np.sin(ang)
        self.body_to_space = np.transpose(self.space_to_body)
