"""
:module: AuCub8
:platform: Unix, Windows
:synopsis: Defines the class for AuCub8

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu>, March2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - added get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division

import numpy as np
from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs


class AuCub8(NanoAbs):
    """Nanocrystal of 8 gold particles, suitable for unit tests
    """

    def __init__(self, ff):
        """

        :param ff: the name of the forcefield to be used to construct this core
        """

        super(AuCub8, self).__init__(ff, 9, 'AuCub8')

        # number of particles in entire system by type (including interior, unsimulated particles)
        self.num_total_particles = np.array([8, 6], dtype=int)

        self.position = np.array([[0, 0, 0], [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1],
                                  [-1, 1, 1], [1, -1, 1], [1, 1, -1], [1, 1, 1]])

        #  this is hardcoded, but it shouldn't be
        self.mass[1:1+8] = self.ff_reader.get_molecular_weight('Au')
        self.mass[0] = np.sum(self.mass[1:])

        # types
        self.types = ['_AuCub8', 'Au']
        self.typeid = ['_AuCub8'] + ['Au']*8
        self.graft_sites = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0], [-1, 0, 0], [0, 0, -1]])
        self.graft_num = 6

        self.moment_inertia[0] = np.array([3151.46512, 3151.46512, 3151.46512])

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[2] - self.position[1]  # points along the z-axis
