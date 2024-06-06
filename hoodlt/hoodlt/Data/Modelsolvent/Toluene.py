"""
:module: Toluene
:platform: Unix, Windows
:synopsis: Class for Toluene solvent

..moduleauthor: Tommy Waltmann <tomwalt@iastate.edu> June 2018
..history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Removed unecessary functions
..              Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                - edited get_vector() method
..                - Updated class to work with BasicSystemEntity
"""

import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class Toluene(SolventAbs):
    """
    Class that implements a Toluene solvent
    """
    def __init__(self, ff):
        """

        :param ff: the name of the forcefield used to construct this object
        """

        super(Toluene, self).__init__(1, ff, 8, 'Toluene')  # always 1 repeat

        # masses, in amu
        self.mass[1] = self.ff_reader.get_molecular_weight('C')
        self.mass[2:7] = self.ff_reader.get_molecular_weight('CH')
        self.mass[7] = self.ff_reader.get_molecular_weight('CH3')
        self.mass[0] = np.sum(self.mass[1:])

        # set the positions
        self._set_positions()

        # types, typeids, etc
        self.types = ['_Toluene', 'C', 'CH', 'CH3']
        self.typeid = ['_Toluene'] + ['C'] + ['CH']*5 + ['CH3']

        # rigid body stuff
        self.moment_inertia[0] = np.diag(self.moment_of_inertia())
        self.body = np.zeros(self.num_particles)

    def _set_positions(self):
        """
        Sets all the positions for the Toluene molecule

        :return: An 8x3 numpy array which contains the positions of the particles on the Toluene
        """

        self.position = np.array([[0, 0, 0], # this position is a dummy, will reset once we shift center of mass to origin
                              [0, 0, 0],
                              [1.4 * np.sin(np.pi / 3), -1.4 * np.cos(np.pi / 3), 0],
                              [1.4 * np.sin(np.pi / 3), -1.4 * (1 + np.cos(np.pi / 3)), 0],
                              [0, -1.4 * (1 + 2 * np.cos(np.pi / 3)), 0],
                              [-1.4 * np.sin(np.pi / 3), -1.4 * (1 + np.cos(np.pi / 3)), 0],
                              [-1.4 * np.sin(np.pi / 3), -1.4 * np.cos(np.pi / 3), 0],
                              [0, 1.54, 0]])

        self.shift(-1 * self._center_of_mass_ignoring_center())  # now the center of mass is at the origin
        self.position[0] = np.array([0.0, 0.0, 0.0])

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]

    def _center_of_mass_ignoring_center(self):
        """
        Calculate the center of mass, ignoring the fictional rigid center.

        :return: a numpy array containing the position of the center of mass
        """

        r_cm = np.zeros(3)
        for i in range(len(self.mass) - 1):
            for j in range(3):
                r_cm[j] += self.mass[i + 1]*self.position[i + 1, j]

        return r_cm / np.sum(self.mass[1:])

