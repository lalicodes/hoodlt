"""
:module: AnalyzeCore
:Platform: Windows, Unix
:synopsis: Class which contains methods to do calculations with individual core objects

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Made to inherit from abstract class
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions
from hoodlt.Analysis.Analyze.AnalyzeBasicEntity import AnalyzeBasicEntity


class AnalyzeCore(AnalyzeBasicEntity):
    """
    Class which contains methods to do calculations with individual core objects
    """

    def __init__(self, core, l_b):
        """

        :param core: NanoAbs object taken from a trajectory
        :param l_b: matrix defining the box
        """

        super(AnalyzeCore, self).__init__(core)

        self.box_dim = BoxFunctions(l_b)

    def radius_area(self):
        """
        Returns surface area of the core, assuming that it is a sphere with radius equal to the average radius of the
        core

        :return: surface area of the core, while making some assumptions
        :rtype: float
        """

        return self.radius_avg() ** 2 * 4 * np.pi

    def radius_avg(self):
        """
        Returns the average distance between the core center and the grafting positions of the core

        :return: average radius of the core
        :rtype: float
        """

        return np.average([self.box_dim.compute_distances(i, self.entity.position[0])[0] for ind, i
                           in enumerate(self.entity.graft_sites)])

    def max_radius(self):
        """
        Returns the largest distance between the core center and any point on the core/grafting site

        :return: value of maximum radius
        :rtype: float
        """

        mat = [self.box_dim.compute_distances(i, self.entity.position[0])[0] for i in self.entity.position]

        return np.max(mat)

    def min_radius(self):
        """
        Returns the smallest distance between the core center and a grafting site

        :return: the minimum radius of the grafters
        """

        mat = [self.box_dim.compute_distances(i, self.entity.position[0])[0] for i in self.entity.graft_sites]

        return np.min(mat)

    def particle_1_spherical_coordinates(self):
        """
        Gets the direction of the first particle on the core's shell in the spherical coordinates theta and phi

        :return: tuple (theta, phi)
        """

        x, y, z = self.entity.position[1] - self.entity.position[0]

        xsqplusysq = x ** 2 + y ** 2
        # r = np.sqrt(XsqPlusYsq + z ** 2)  # r
        phi = np.pi / 2 + np.atan2(z, np.sqrt(xsqplusysq))  # phi
        theta = np.atan2(y, x)  # theta
        return theta, phi
