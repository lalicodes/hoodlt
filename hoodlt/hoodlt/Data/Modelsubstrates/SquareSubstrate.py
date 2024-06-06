"""
:module: SubstrateAbs
:platform: Unix, Windows
:synopsis: Test system to use for a substrate

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov> June 2022
..                  - changed the organization of the class
..                  - Made it consistent with HOOMD v3
"""

import numpy as np
from hoodlt.Data.Modelsubstrates.SubstrateAbs import SubstrateAbs


class SquareSubstrate(SubstrateAbs):
    """
    A substrate which is built in the x-y plane and is centered at the origin, in a square lattice style.
    """

    def __init__(self, forcefield, grid_dim, spacing, particle_type='Au'):
        """

        :param forcefield: the forcefield to built this substrate with
        :param grid_dim: the number of particles of the grid in the x and y directions, input is of the form [x, y]. Must be an odd number so the rigid center is at the origin
        :param spacing: the distance between adjacent particles on the substrate
        :param particle_type: particle type to use for the particles in the grid
        """

        if (grid_dim[0] + 1) % 2 == 0 or (grid_dim[1] + 1) % 2 == 0:
            raise ValueError("Both grid dimensions must be an even number")

        self.num_grid_points_x = grid_dim[0]
        self.num_grid_points_y = grid_dim[1]

        name = particle_type + str(self.num_grid_points_x) + "-" + str(self.num_grid_points_y)

        super(SquareSubstrate, self).__init__(forcefield, 1 + self.num_grid_points_x*self.num_grid_points_y, name)

        self.types = ['_'+particle_type, particle_type]
        self.typeid = [self.types[0]] + [self.types[1]]*(self.num_particles - 1)

        grid_pos_x = self._generate_positions(spacing, self.num_grid_points_x)
        grid_pos_y = self._generate_positions(spacing, self.num_grid_points_y)

        pos = [[0.0, 0.0, 0.0]]
        for pos_x in grid_pos_x:
            for pos_y in grid_pos_y:
                pos.append([pos_x, pos_y, 0.0])

        # set the positions
        self.position = np.array(pos)

        # set the masses
        self.mass[1:self.num_particles] = self.ff_reader.get_molecular_weight(particle_type)
        self.mass[0] = self.ff_reader.get_molecular_weight(particle_type) * self.num_particles

        self.moment_inertia[0] = np.diag(self.moment_of_inertia())
        self.periodic_box = np.array([grid_dim[0], grid_dim[1], np.inf, 0.0, 0.0, 0.0])

    @staticmethod
    def _generate_positions(spacing, num_points):
        """
        Generate a list of points of distance spacing between each other num_points long

        :param spacing: the spacing between the points
        :param num_points: the number of points to generate
        :return: a list of points each separated from each other by spacing
        """
        num = num_points//2
        return [(0.5+ind)*spacing for ind in range(-num, num)]
