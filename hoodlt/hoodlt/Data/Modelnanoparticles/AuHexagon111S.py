"""
:module: Au111S
:platform: Unix, Windows
:synopsis: Defines the class for AuHexagon111S

.. moduleauthor:: Jacob Austin <jaustin2@iastate.edu>, June 2019
.. history::
..                Jacob Austin <jaustin2@iastate.edu> April 2020
                    - added special coulomb and lj functionality
..
..
"""

from __future__ import division

import numpy as np
from hoodlt.Data.Modelnanoparticles.CoreWithPositionsInFile import CoreWithPositionsInFile


class AuHexagon111S(CoreWithPositionsInFile):
    """
    Defines the AuHexagon111S core object
    """
    num_core_particles = 91
    def __init__(self, ff, core_types=['Au'], scale=1, name_rigid_center=None):

        """
        :param ff: the name of the forcefield to be used to construct this object
        :param core_types: list of types of the core particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting
        sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the
        typical '_' prefix
        """

        super(AuHexagon111S, self).__init__(ff, self.num_core_particles, "AuHexagon111S.txt", core_types, scale, name_rigid_center)

        self.moment_inertia[0] = np.array([9008.6042, 9043.2008, 18051.8050])


    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[46] - self.position[0]  # points along the z-axis



