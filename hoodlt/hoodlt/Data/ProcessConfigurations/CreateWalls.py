"""
:module: CreateWalls
:platform: Unix, Windows
:synopsis: Defines walls to be used by squeeze

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> June 2022
.. history:
..                Jonas Hallstrom <jonasleo@hotmail.com> June 2022
..                  - fixed the sign of the normal vectors of the planar walls
"""


import numpy as np
import numpy.linalg as la
import hoomd


class CreateWalls:
    def __init__(self, name):
        """Takes a list of FunctionalizedParticle objects and squeezes each one into the desired radius and shape

        :param name: wall name
        """

        # use lower case names
        name.lower()
        planar = ['plane', 'slab', 'square', 'cube']
        self.wall_names = ['sphere', 'cylinder'] + planar

        if name not in self.wall_names:
            print('Available walls are ', self.wall_names)
            print(name + ' is not on the list of available names')
            raise ValueError('incorrect wall name')

        self.wall_name = name
        # define the norm for sphere and cylinder
        self.order = None
        self.num = 3
        # for cylinder, only the x,y coordinates need to be accounted for
        if self.wall_name == 'cylinder':
            self.num = 2

        self.planar_dict = {'plane': 1, 'slab': 2, 'square': 4, 'cube': 6}
        # normal to planes
        self.normals = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0,  1, 0), (0, 0, -1), (0, 0,  1)]
        if self.wall_name in planar:
            self.index_planar = self.planar_dict[self.wall_name]
            self.wall_name = 'planar'
            self.num = self.index_planar//2
            self.order = np.inf

    def make_wall(self, size):
        """
        returns the wall with the given size
        """
        return getattr(self, self.wall_name)(size)

    def max_radius(self, mat):
        """
        returns the largest distance to the wall
        """

        return np.max(getattr(self, 'def_norm')(mat, self.num, self.order))

    @staticmethod
    def sphere(size):
        """
        sphere wall

        :param size: sphere radius
        :return: sphere wall
        """

        return [hoomd.wall.Sphere(radius=size)]

    @staticmethod
    def cylinder(size):
        """
        sphere wall

        :param size: cylinder radius
        :return: cylinder wall
        """

        return [hoomd.wall.Cylinder(radius=size)]

    def planar(self, size):
        """
        sphere wall

        :param size: plane distance
        """

        lst = []
        for ind in range(self.index_planar):
            # this defines a point in the plane
            orgn = tuple(-size*np.array(self.normals[ind]))
            lst.append(hoomd.wall.Plane(origin=orgn, normal=self.normals[ind]))

        return lst

    @staticmethod
    def def_norm(vec, num, order):
        """
        :param vec: numpy array
        :param num:
        :param order: the type of norm to be used
        :return norm
        """

        # plane is special
        if num == 0:
            return vec[:, 0]

        return la.norm(vec[:, :num], axis=1, ord=order)
