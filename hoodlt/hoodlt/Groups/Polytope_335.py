"""
:module: Polytope_335
:platform: Unix, Windows, OS
:synopsis: Defines the Polytope 335

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  May 2017
"""

from __future__ import division
import numpy as np
import hoodlt.Groups.Quaternion as Quat
import hoodlt.Groups.Yico_Group as YcG


class Polytope335(object):
    """
    Defines the polytope :math:`\\{ 3,3,5 \\}` in four dimensional space, using quaternionic representation
    The polytope includes the point (1, 0, 0, 0) as one of the vertices
    
    all vertices are inscribed in a hypersphere of radius 1 and edge separation :math:`1/\\tau`
    
    :math:`\\tau` is the golden ratio.
    
    """

    def __init__(self, radius_s3=1):
        """The Constructor
        
        :param radius_s3: Radius of the hypersphere
        """

        # radius
        self.radius = radius_s3

        # icosahedral group
        self.ico = YcG.YGroup()

        # vertices
        self.num_vertices = 120
        v_poly = np.zeros([self.num_vertices, 4])
        self.vertices = Quat.Quaternion(v_poly)

        # reference point
        self.tau = 0.5*(1+np.sqrt(5))
        vref = self.radius*np.array([1.0, 0, 0.0, 0.0])
        qref = Quat.Quaternion(vref)

        num_cum = 0
        for enum_class, g_class in enumerate(self.ico.classes):
            for ind in range(g_class):
                qy = self.ico.group_class(enum_class, ind)
                qr = qy.inverse() * qref
                self.vertices.quat[num_cum] = qr.quat
                num_cum += 1
