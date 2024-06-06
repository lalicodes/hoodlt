"""
:module: Water
:platform: Unix, Windows
:synopsis: Implements a class defining water solvents

.. moduleauthor:: Elizabeth Macias
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class SpceWater(SolventAbs):
    """
    Defines SPC/E model H2O as a rigid body
    """

    def __init__(self, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        super(SpceWater, self).__init__(1, ff, 3, 'SpceWater')


        # types (_Water, O, H)
        self.types = ['OW', 'H']
        self.typeid = ['OW']*1 + ['H']*2

        self.diameter = np.zeros(self.num_particles) + 0.08

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('OW') # amu
        self.mass[1] = self.ff_reader.get_molecular_weight('H') # amu
        self.mass[2] = self.ff_reader.get_molecular_weight('H') # amu


        # charges
        self.charge[0] = self.ff_reader.get_charge('OW')
        self.charge[1] = self.ff_reader.get_charge('H')
        self.charge[2] = self.ff_reader.get_charge('H')

        # defining constituent particle positions
        theta_spce = self.ff_reader.get_angle_t0('H-OW-H') # H-O-H angle in radians
        roh = self.ff_reader.get_bond_r0('OW-H') # distance between H and O.

        # particle positions
        self.position[0] = [0, 0, 0] # O, H, and H positions
        self.position[1] = [roh*np.sin(theta_spce/2), roh*np.cos(theta_spce/2), 0]
        self.position[2] = [roh*np.sin(-theta_spce/2), roh*np.cos(theta_spce/2), 0]

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[-1] - self.position[0]

    def get_constraints(self):
        """
        override of method in BasicSystemEntity
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system

        :return: a 2x2 numpy array
        """

        c_list = [[0, 1], [0, 2]]

        return np.array(c_list, dtype=int)

    def get_constraint_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the bond types in the system
        """

        c_types = ['OW-H'] * 2

        dist_c = []
        for ctyp in c_types:
            dist_c.append(self.ff_reader.get_bond_r0(ctyp))

        return dist_c

    def get_angles(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (3 particles per angle) involved in the angles in the system

        :return: a 1x3 numpy array
        """

        return np.array([[1, 0, 2]])


    def get_angle_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the angle types in the system
        """

        return ['H-OW-H']
