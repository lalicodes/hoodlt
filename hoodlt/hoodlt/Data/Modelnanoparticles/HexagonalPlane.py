"""
:module: AuHexagon71S
:platform: Unix, Windows
:synopsis: Defines the class for AuHexagon71S

.. moduleauthor:: Jacob Austin <jaustin2@iastate.edu>, June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - made the class inherit from the new abstract class, and renamed it
..
..
"""

from __future__ import division

import numpy as np
from hoodlt.Data.Modelnanoparticles.CoreWithPositionsInFile import CoreWithPositionsInFile


class HexagonalPlane(CoreWithPositionsInFile):
    """
    Defines the HexagonPlane core object
    """

    num_core_particles = 61

    def __init__(self, ff, core_types=['Au'], scale=1, name_rigid_center=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param core_types: list of types of the core particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting
        sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the
        typical '_' prefix
        """

        super(HexagonalPlane, self).__init__(ff, self.num_core_particles, "HexagonalPlane.txt", core_types, scale, name_rigid_center)

        self.moment_inertia[0] = np.array([4004.431925, 4019.3763, 8011.2407])

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[31] - self.position[0]  # points along the z-axis
