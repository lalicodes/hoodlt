"""
:module: AnalyzeSolvent
:Platform: Windows, Unix
:synopsis: Class which contains methods to do calculations with individual solvent objects

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Made to inherit from abstract class
"""

import numpy as np
from hoodlt.Analysis.Analyze.AnalyzeBasicEntity import AnalyzeBasicEntity


class AnalyzeSolvent(AnalyzeBasicEntity):
    """
    Class which contains methods to do calculations with individual solvent objects
    """

    def __init__(self, solv):
        """

        :param solv: SolventAbs object
        """

        super(AnalyzeSolvent, self).__init__(solv)

    def radius_of_gyration(self):
        """
        Calculate the radius of gyration for a solvent molecule

        :return:
        """

        n = self.entity.position.shape[0]
        rog = []
        com = self.entity.center_of_mass()

        for i in range(n):
            r = np.subtract(self.entity.position[i], com)**2
            rog.append(r)

        rad_gyr = np.sqrt((1/n)*np.sum(rog))

        return rad_gyr

    def monomer_radius_of_gyration(self, num_monomers, len_monomers):
        """
        Calculate the monomer radius of gyration for a solvent molecule.
        Note that self.num_particles = num_monomers*len_monomers

        :param num_monomers: number of monomers on the solvent
        :param len_monomers: length (number of particles) of a monomer on the solvent
        :return: the monomer radius of gyration of the solvent
        """

        n = self.entity.position.shape[0]/len_monomers

        print(len_monomers)
        print(num_monomers)

        rog = []
        com = self.entity.center_of_mass()
        for i in range(int(n)):
            r = np.subtract(self.entity.position[i*len_monomers], com)**2
            rog.append(r)

        rad_gyr = np.sqrt((1/n)*np.sum(rog))

        return rad_gyr
