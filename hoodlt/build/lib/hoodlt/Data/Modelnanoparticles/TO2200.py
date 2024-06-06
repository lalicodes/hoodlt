"""
:module: TO2200
:platform: Unix, Windows
:synopsis: Defines the class for TO2200 (AuShell2200S)

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu>, October 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - class can have any core types
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - added get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division

from hoodlt.Data.Modelnanoparticles.CoreWithPositionsInFile import CoreWithPositionsInFile


class TO2200(CoreWithPositionsInFile):
    """
    Defines a truncated octahedron with 2200 core particles on the shell, and 750 grafting sites cut from the 1-1-1 and
    0-0-1 planes of an FCC lattice
    """

    num_core_particles = 2200

    def __init__(self, ff, core_types=['Au'], scale=1, name_rigid_center=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param core_types: list of types of the core particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the typical '_' prefix
        """

        super(TO2200, self).__init__(ff, self.num_core_particles, "TO2200.txt", core_types, scale, name_rigid_center)

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]  # points toward the first particle on the core
