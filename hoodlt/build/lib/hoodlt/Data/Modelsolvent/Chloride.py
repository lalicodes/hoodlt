"""
:module: Water
:platform: Unix, Windows
:synopsis: Implements a class defining water solvents

.. moduleauthor:: Elizabeth Macias
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class Chloride(SolventAbs):
    """
    Defines Cl
    """

    def __init__(self, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        super(Chloride, self).__init__(1, ff, 1, 'Chloride')


        # types ('Cl')
        self.types = ['Cl']
        self.typeid = ['Cl']*1

        self.diameter = np.zeros(self.num_particles) + 0.08

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('Cl') # amu

        # charges
        self.charge[0] = self.ff_reader.get_charge('Cl')

        # particle positions
        self.position = np.zeros((1, 3))
