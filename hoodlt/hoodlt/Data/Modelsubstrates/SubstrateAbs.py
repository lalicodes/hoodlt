"""
:module: SubstrateAbs
:platform: Unix, Windows
:synopsis: Defines the abstract classes used model substrates

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov> June 2022
..                  - added member periodic box
..                  - Made it consistent with HOOMD v3

"""

import numpy as np
from hoodlt.Data.Modelconfigurations.BasicSystemEntity import BasicSystemEntity


class SubstrateAbs(BasicSystemEntity):
    """
    Abstract class to model substrates
    """

    def __init__(self, forcefield, num_particles, name):
        """

        :param forcefield: the forcefield used to build the substrate
        :param num_particles: the number of particles in the substrate
        :param name: the name of the substrate
        """

        super(SubstrateAbs, self).__init__(forcefield, num_particles, name)

        self.body = np.zeros(self.num_particles)  # all substrates should be one single rigid body
        # sets the size of the periodic box (to be defined by derived classes)
        self.periodic_box = None

    def get_vector(self):
        """
        Gets a vector normal to the plane of the substrate

        :return: a vector normal to the plane of the substrate
        :rtype: numpy array, [x, y, z]
        """

        # get 2 vectors in the plane of the substrate
        plane_vec_1 = self.position[1] - self.position[0]
        plane_vec_2 = self.position[2] - self.position[0]

        return np.cross(plane_vec_1, plane_vec_2)  # cross product will necessarily be normal to the substrate
