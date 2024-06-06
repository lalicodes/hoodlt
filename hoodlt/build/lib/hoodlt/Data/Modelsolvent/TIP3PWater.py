"""
:module: TIP3PWater
:platform: Unix. Windows
:synopsis: Implements a class defining TIP3P WaterMolecule

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> August 2019
.. history:
"""

from __future__ import division
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs
import numpy as np


class TIP3PWater(SolventAbs):
    """
    Defines solvent  TIP3P Water (:math:`\\mbox{H}_2\\mbox{O}`)
    """

    def __init__(self, ff):
        """
        create the TIP3P Water molecule
        """

        super(TIP3PWater, self).__init__(1, ff, 4, 'TIP3PWater')  # always only 1 repeat

        self.mass[1] = self.ff_reader.get_molecular_weight('OH2')
        self.mass[2:4] = self.ff_reader.get_molecular_weight('H')
        self.mass[0] = np.sum(self.mass[1:])

        self.charge[1] = self.ff_reader.get_charge('OH2')
        self.charge[2:4] = self.ff_reader.get_charge('H')

        # types
        self.types = ['_TIP3PWater', 'OH2', 'H']
        self.typeid = ['_TIP3PWater'] + ['OH2'] + ['H'] * 2

        # particle positions, in Angstroms
        self._set_positions()

        self.moment_inertia[0] = np.diag(self.moment_of_inertia())
        self.body = np.zeros(self.num_particles)

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]  # points from O to the rigid center

    def _set_positions(self):
        """
        called by initializer to get positions

        :return: position array for the water molecule
        """
        self.position = np.zeros((4, 3))

        r_oh = self.ff_reader.get_bond_r0('OH2-H')

        theta_oho = self.ff_reader.get_angle_t0('H-OH2-H')

        self.position[2] = [r_oh*np.cos(theta_oho/2), r_oh*np.sin(theta_oho/2), 0]
        self.position[3] = [r_oh*np.cos(theta_oho/2), -r_oh*np.sin(theta_oho/2), 0]

        self.shift(-1 * self._center_of_mass_ignoring_center())  # now the center of mass is at the origin
        self.position[0] = np.array([0.0, 0.0, 0.0])

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

    def get_bonds(self):
        """
        override of method in BasicSystemEntity
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system

        :return: a 2x2 numpy array
        """

        return np.array([[1, 2], [1, 3]])

    def get_angles(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (3 particles per angle) involved in the angles in the system

        :return: a 1x3 numpy array
        """

        return np.array([[2, 1, 3]])

    def get_bond_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the bond types in the system
        """

        return ['OH2-H'] * 2

    def get_angle_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the angle types in the system
        """

        return ['H-OH2-H']