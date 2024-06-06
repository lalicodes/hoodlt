"""
:module: NanoAbs
:platform: Unix, Windows
:synopsis: Defines the abstract classes used to store nanoparticles

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - updated grafting algorithm and removed some code we don't need anymore
..                  - orienting and de-orienting cores now updates the grafting sites as well
..                  - cores now have a list of graft sites separate from their positions list
..                  - moved many methods to analysis
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - Documentation
..                  - added get_vector for use in aligning cores
..                  - removed references to body_count variable when dumping
..                  - removed the '_' in front of many of the names of the methods
..                  - core_name is now set in basic entity class
..                Tommy Waltmann <tomwalt@iastate.edu> May 2019
..                  - made class inherit from BasicSystemEntity, moved many methods to that class
..                  - removed uncessary values stored by the constructor
..                  - rewrote many methods to be simpler, easier to understand
..                Alex Travesset <trvsst@ameslab.gov> June 2022
..                  - included rotations consistent with v3 (eliminated old orient)
"""

import numpy as np
import numpy.linalg as la
from scipy.spatial import ConvexHull
from hoodlt.Data.Modelconfigurations.BasicSystemEntity import BasicSystemEntity
from hoodlt.Groups.Quat import Quat


class NanoAbs(BasicSystemEntity):
    """
    Defines general nanoparticle cores
    """

    def __init__(self, forcefield, num_particles, name):
        """

        :param forcefield: the name of the forcefield to be used to construct this object
        :param num_particles: the number of particles (including rigid center) on this core
        :param name: the name of the core
        """

        super(NanoAbs, self).__init__(forcefield, num_particles, name)

        # grafter stuff
        self.graft_sites = np.array([])
        self.graft_num = None

        self.body = np.zeros(self.num_particles)  # each core will be one rigid body
        self.params_outside_forcefield = {}

    def volume(self):
        """Computes the nanoparticle volume

        :return:
        :rtype: float
        """

        raise ValueError("volume() method in NanoAbs is not implemented")

    def area(self):
        """returns the area using formula

        :return: value of the area
        :rtype: float
        """

        raise ValueError("area() method in NanoAbs is not implemented")

    def rotate(self, quat):
        """
        rotates the core as well as its grafting sites
        (this is used in preparing simulations)

        :param quat: quaternion as a 4-tuple
        :return:
        """
        qt = Quat()
        quatp = self.orientation[0]

        self.orientation[0] = tuple(qt.multiply(quat, quatp))

    def rotate_actual(self, quat):
        """
        rotates the core as well as its grafting sites and implements the rotation
        (this is used to draw actual rotated nanoparticles)

        :param quat: quaternion as a 4-tuple
        :return:
        """
        qt = Quat()
        relative_center = self.position[0]
        for ind in range(len(self.position)):
            self.position[ind] = relative_center + qt.rotation(quat, self.position[ind]-relative_center)

    def align_core(self, vector):
        """
        Aligns the core's instrinsic vector with the input vector, using the rigid center as the origin for the rotation

        :param vector: the vector to align the core's intrinsic vector with
        :return: None
        """

        self.align(self.position[0], vector)

    def shift(self, vector):
        """
        Shift the core by the given vector

        :param vector: the vector to shift the core by
        :return:
        """

        super(NanoAbs, self).shift(vector)

        # shift the graft sites as well
        self.graft_sites = self.graft_sites[:] + vector

    def align_core_using_vectors(self, original_vector, vector_to_rotate_to):
        """
        Rotates the core in the exact same way that original_vector must be rotated so that it points in the same
        direction as vector_to_rotate_to
        :param original_vector: vector that will be used as the start reference for the rotation of the entity
        :param vector_to_rotate_to: vector that will serve as end reference for the rotation of the entity
        :return: None
        """

        self.align_using_vectors(self.position[0], original_vector, vector_to_rotate_to)

    def _project_on_surface(self, point):
        """
        
        :param point: som point in space
        :return: position on core
        """

        unit = point / la.norm(point)

        # get the three points on the core which are the closest unit vector to the input point
        shell_positions = self._get_points_on_core_surface()
        vectors = np.array([np.dot(unit, pos / la.norm(pos)) for pos in shell_positions])
        arg = np.argsort(vectors)[::-1]
        closest_0 = self.position[arg[0]]
        closest_1 = self.position[arg[1]]
        closest_2 = self.position[arg[2]]

        # scale the input point to be the length of the closest point
        new_point = la.norm(closest_0) * unit

        return new_point

    def _get_points_on_core_surface(self):
        """
        Gets some points which lie on the surface of the core, by computing its convex hull

        :return: A numpy array of the points on the core which constitute its surface
        """

        shell_indexes = ConvexHull(self.position).vertices
        shell_positions = self.position[shell_indexes]

        return shell_positions

    @staticmethod
    def points_on_unit_sphere(n):
        """
        Calculates n points distributed on unit sphere

        :param n: number of points to distribute
        :return: positions of points on unit sphere
        """

        pts = []
        inc = np.pi * (3 - np.sqrt(5))
        off = 2 / float(n)
        for k in range(n):
            y = k * off - 1 + (off/2)
            r = np.sqrt(1 - y*y)
            phi = k * inc
            pts.append([np.cos(phi)*r, y, np.sin(phi)*r])
        return np.array(pts)

    def add_grafters(self):
        """
        Creates grafting sites on the core

        :return:
        """

        pnts = self.points_on_unit_sphere(self.graft_num)
        self.graft_sites = np.array([self._project_on_surface(y) for y in pnts])

    def get_name(self):
        """
        name

        :returns: name
        """

        return self.name

    def parameters_outside_force_field(self):
        """
        parameters defined in the core that do no follow from the force field

        """

        return self.params_outside_forcefield
