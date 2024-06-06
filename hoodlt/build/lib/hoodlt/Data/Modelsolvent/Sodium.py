"""
:module: Water
:platform: Unix, Windows
:synopsis: Implements a class defining water solvents

.. moduleauthor:: Elizabeth Macias
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class Sodium(SolventAbs):
    """
    Defines Na
    """

    def __init__(self, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        super(Sodium, self).__init__(1, ff, 1, 'Sodium')


        # types ('Na')
        self.types = ['Na']
        self.typeid = ['Na']*1

        self.diameter = np.zeros(self.num_particles) + 0.08

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('Na')

        # charges
        self.charge[0] = self.ff_reader.get_charge('Na')

        # particle positions
        self.position = np.zeros((1, 3))
