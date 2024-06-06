"""
:module: PolymerLigand
:platform: Unix, Windows
:synopsis: Implements the polymer version of the abstract ligand class

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> May 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs

class PolymerLigand(LigandAbs):
    """Defines the Polymer Ligand
    """

    def __init__(self, repeats, ff, helix=None):
        """The constructor

        :param repeats: number of repeatable polymer units in chain
        :param helix: optional length to pass for helix creation
        """

        super(PolymerLigand, self).__init__(repeats, ff, repeats+1, 'Polymer')

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('Gr')
        self.mass[1:repeats + 1] = self.ff_reader.get_molecular_weight('O')

        # types
        self.types = ['Gr', 'O']
        self.typeid = ['Gr']*1 + ['O']*repeats

        self.bond_length = 1

        if helix is None:
            self.helix_length = self.num_particles * self.bond_length
        else:
            self.helix_length = helix

        # positions
        self.position = self.build_chain()

    def build_chain(self):
        """builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """
        chain = np.zeros((np.sum(self.num_particles), 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polybond')

        current = np.array([-1.0, 0.0, 0.0])

        for i in range(0, self.num_particles):

            pt = np.array([1.0, 0.0, 0.0])*float(poly_bond_length)

            current += pt

            chain[i] = current

        sub_pt = np.array([0, 0, 0])
        for i in range(self.num_particles):
            chain[i] = np.subtract(chain[i], sub_pt)

        return chain

    def get_bond_types(self):
        """utilized to dump_gsd

        :return: list of all the bond types in the system
        """
        return ['polybond'] * self.repeats

    def get_helix_point(self, length, radius=0.9, pitch=1, turns=2):
        """Output point position along helix

        :param length: length along helix
        :param radius: radius of helix coil
        :param pitch: pitch of helix coil
        :param turns: number of turns in helix
        :return: position of helix point
        """

        z = radius * (np.cos(length * turns) - 1)
        y = radius * np.sin(length * turns)
        x = pitch * length

        return np.array([x, y, z])

    def get_helical_radius(self, arc_length, dist):
        """Get helix radius from given arc_length and helix coil height

        :param arc_length: length of arc constituting helix
        :param dist: height of helix coil
        :return: helix radius
        """

        return np.sqrt(arc_length**2 - dist**2) / (2 * np.pi)

    def get_helical_length(self, dist, radius):
        """Get helical arc length from given helix coil height and radius

        :param dist: height of helix coil
        :param radius: radius of helix
        :return: helical arc length
        """

        return np.sqrt(dist**2 + (2 * np.pi * radius)**2)

    def get_angles(self):
        return []

    def get_dihedrals(self):
        return []
