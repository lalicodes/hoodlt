"""
:module: TO201
:platform: Unix, Windows
:synopsis: Defines the class for TO201 (Au201S)

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu>, March2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - inherits from the new abstract class
..                  - renamed associated text file with positions and grafting sites
..                  - class is now generalized and can handle multiple particles types on the core, as well as scaling
..                    the core positions and grafting sites
..                  - Au201S is now TO201
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - added get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division

from hoodlt.Data.Modelnanoparticles.CoreWithPositionsInFile import CoreWithPositionsInFile


class TO201(CoreWithPositionsInFile):
    """
    Defines a truncated octahedron with 122 core particles on the shell, and 80 grafting sites cut from the 1-1-1 and
    0-0-1 planes of an FCC lattice
    """

    num_core_particles = 122

    def __init__(self, ff, core_types=['Au'], scale=1, name_rigid_center=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param core_types: list of types of the core particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the typical '_' prefix
        """

        super(TO201, self).__init__(ff, self.num_core_particles, "TO201.txt", core_types, scale, name_rigid_center, 201)

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[71] - self.position[0]  # points along the z-axis
