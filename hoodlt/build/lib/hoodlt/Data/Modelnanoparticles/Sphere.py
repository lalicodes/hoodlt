"""
:module: Sphere
:platform: Unix, Windows
:synopsis: Defines the class for Spherical particles

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu>, April 2017
.. history:
                Tommy Waltmann <tomwalt@iastate.edu> June 2019
                  - added get_vector() method
                  - reworked class so it fits with BasicSystemEntity
                  - name is now set in basic entity class
                Xun Zha <xzha@iastate.edu> June 2021
                    - Reading of the file consistent with the units in the forcefield (adding Angstrom to Construction)
"""


from __future__ import division

import numpy as np
from scipy import spatial as sp

from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs


class Sphere(NanoAbs):
    """Defines a general spherical nanoparticle
    """
    def __init__(self, sphere_info, ff):
        """
        Sphere core

        :param sphere_info: SphereInfo object, containing: (core radius, grafting radius, number of atoms in core,
         number of grafting sites)
         :param ff: the name of the forcefield to be used to construct this core
        """

        # sphere information
        rc = sphere_info.core_radius
        rg = sphere_info.graft_radius
        nc = sphere_info.core_num
        ng = sphere_info.graft_num

        # central particle type
        core_name = 'Sphere-' + str(int(rg)) + '-Gr-' + str(ng)

        super(Sphere, self).__init__(ff, nc+1, core_name)

        self.mass[1:1+nc] = self.ff_reader.get_molecular_weight('Cr')
        self.mass[0] = np.sum(self.mass[1:])

        # types
        self.types = ['_'+self.name, 'Cr']
        self.typeid = ['_'+self.name] + ['Cr']*nc
        self.graft_num = ng

        # get the units of the forcefield
        angstrom_to_construction = self.ff_reader.get_units().angstrom_to_construction

        # particle positions
        self.pos_core = np.array(self.points_on_unit_sphere(nc))*rc * angstrom_to_construction
        self.graft_sites = np.array(self.points_on_unit_sphere(ng))*rg * angstrom_to_construction
        self.position = np.array([[0.0, 0.0, 0.0]] + self.pos_core.tolist())

        print('CORE: MEAN / MIN / MAX')
        self.sphere_nn(self.pos_core)
        print('GRAFT: MEAN / MIN / MAX')
        self.sphere_nn(self.graft_sites)

        # moment of inertia in rest frame
        # Taken directly from NanoAbs class because each instance of sphere will be different
        # need a low tolerance because sphere is not very diagonal
        self.moment_inertia[0] = np.diag(self.moment_of_inertia())

    def sphere_nn(self, pts):
        """calculates NN distance of particles on Sphere

        :param pts: points to consider for NN operation
        :return: average NN distance, minimum NN distance, maximum NN distance
        """

        nn = sp.KDTree(pts)
        nn1 = nn.query(pts, k=2)[0][:, 1]
        mean_dist, min_dist, max_dist = np.mean(nn1), np.amin(nn1), np.amax(nn1)
        print(mean_dist, '/', min_dist, '/', max_dist)

        return mean_dist, min_dist, max_dist

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]  # points toward the position of the first core particle
