"""
:module: AnalyzeNanoparticle
:Platform: Windows, Unix
:synopsis: Class which contains methods to do calculations with individual nanoparticle objects

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Corrected an import statement
"""

import numpy as np
import numpy.linalg as la
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions
from hoodlt.Analysis.Analyze.AnalyzeLigand import AnalyzeLigand


class AnalyzeNanoparticle:
    """
    This class contains methods to do calculations with individual nanoparticle objects
    """

    def __init__(self, nc, l_b):
        """

        :param nc: FunctionalizedParticle object taken from a specific frame of a trajectory object
        :param l_b: matrix defining the box
        """

        self.nc = nc
        self.box_dim = BoxFunctions(l_b)

    def avg_distance_nth_particle(self, n):
        """
        Computes the average distance between the nth particle on the ligand chain and the core center of this
        nanoparticle

        :param n: nth particle in the ligands aka chain_pos. grafting particle is indexed as 0
        :return: the average distance to the nth particle for all the ligands on the chains
        """

        distance = 0
        for lig in self.nc.ligands:
            if len(lig.position) - 1 >= n:
                distance += AnalyzeLigand(lig).distance_nth_particle(n)
        return distance / len(self.nc.ligands)

    def intermolecular_angles(self):
        """

        :return: a list of all the angles between all the chain vectors in the system
        """

        chain_vecs = [lig.get_vector() for lig in self.nc.ligands]
        inter_angles = []
        for i in range(len(chain_vecs)):
            for x in range(i, len(chain_vecs)):
                g = np.dot(chain_vecs[i], chain_vecs[x]) / (la.norm(chain_vecs[i]) * la.norm(chain_vecs[x]))
                if np.abs(g - 1.0) < 10e-5:
                    g = .999
                inter_angles.append(np.degrees(np.arccos(g)))
        return inter_angles

    def radius_gyration(self):
        """

        :return: radius of gyration of the particle
        """

        rg = np.mean([self.box_dim.compute_distances(lig.position[0], self.nc.core.position[0])[0]**2 for
                        lig in self.nc.ligands])

        return np.sqrt(rg)

    def max_radius(self):
        """
        Determines the maximum distance of any particle on any ligand from the rigid core center

        :param box_dim: optional argument which is used to account for nanoparticles which may be split across the box. If not left to default, enter the [x, y, z] box dimensions
        :return: the largest distance between the center and a ligand particle
        """

        maxes = [0]
        for lig in self.nc.ligands:
            for pos in lig.position:
                dist = self.box_dim.compute_distances(pos, self.nc.core.position[0])[0]
                maxes.append(dist)

        return np.amax(maxes)

    def average_hydrodynamic_radius(self):
        """
        Determines the average distance of the last particle on each ligand from the center of the rigid core.

        :param box_dim: optional argument which is used to account for nanoparticles which may be split across the box. If not left to default, enter the [x, y, z] box dimensions
        :return: the average distance from the center to the last ligand in every chain
        """

        maxes = [0]
        for lig in self.nc.ligands:
            for pos in lig.position:
                dist = self.box_dim.compute_distances(pos, self.nc.core.position[0])[0]
                maxes.append(dist)

        return np.average(maxes)
