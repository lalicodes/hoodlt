"""
:module: AnalyzeLigand
:Platform: Windows, Unix
:synopsis: Class which contains methods to do calculations with individual ligand objects on a nanoparticle

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Made to inherit from abstract class
"""

import numpy as np
import numpy.linalg as la
from hoodlt.Analysis.Analyze.AnalyzeBasicEntity import AnalyzeBasicEntity


class AnalyzeLigand(AnalyzeBasicEntity):
    """
    Class which contains methods to do calculations on specific ligand objects
    """

    def __init__(self, lig):
        """

        :param lig: LigandAbs object, taken from a trajectory
        """

        super(AnalyzeLigand, self).__init__(lig)

    def chain_length(self):
        """
        Computes the length of the chain

        :return: value of the length
        :rtype: float
        """

        return la.norm(self.entity.position[0]-self.entity.position[self.entity.get_num_particles() - 1])

    def distance_nth_particle(self, n):
        """
        Calculates the distance from the nth particle on the chain to the first or 0th particle in the chain

        :param n: nth particle on the chain, an integer
        :return: distance from the 0th to nth particle on the chain
        """

        if n >= self.entity.get_num_particles() or n < 0:
            raise ValueError("The chain does not contain a " + str(n) + "th particle. Indexing starts at 0.")

        return la.norm(np.subtract(self.entity.position[0], self.entity.position[int(n)]))

    def angle_nth_particle(self, n):
        """
        Computes the spherical coordinates between the grafting site and nth particle of the ligand
        :param n: nth particle in the ligands, aka lig_pos. grafting particle is indexed as 0
        :return: the spherical coordinates for the angles between the grafting particle and the nth particle
        in the ligand
        """
        if n >= self.entity.get_num_particles() or n < 0:
            raise ValueError("The chain does not contain a " + str(n) + "th particle. Indexing starts at 0.")

        vector = self.entity.position[int(n)] - self.entity.position[0]

        r = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
        theta = np.arccos(vector[0] / r)
        if vector[0] == 0:
            psi = 0
        else:
            psi = np.arctan(vector[1] / vector[0])

        spherical_coords = np.array([r, theta, psi])
        return spherical_coords
