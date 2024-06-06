"""
:module: TO1289
:platform: Unix, Windows
:synopsis: Defines the class for TO1289 (Au1289S)

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu>, June 2017
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


class TO1289(CoreWithPositionsInFile):
    """
    Defines a truncated octahedron with 482 core particles on the shell, and 258 grafting sites cut from the 1-1-1 and
    0-0-1 planes of an FCC lattice
    """

    num_core_particles = 482

    def __init__(self, ff, core_types=['Au'], scale=1, name_rigid_center=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param core_types: list of types of the core particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the typical '_' prefix
        """

        super(TO1289, self).__init__(ff, self.num_core_particles, "TO1289.txt", core_types, scale, name_rigid_center, 1289)

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[257] - self.position[0]  # points along the z-axis
