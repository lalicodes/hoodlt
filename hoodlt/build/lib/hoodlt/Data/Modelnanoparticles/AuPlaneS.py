"""
:module: AuPlaneS
:platform: Unix, Windows
:synopsis: Defines the class for AuPlaneS

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu>, November 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - added get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from __future__ import division

import numpy as np
import copy as cp
import numpy.linalg as la
from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs


class AuPlaneS(NanoAbs):
    """
    Defines the Au201S nanoparticle
    """

    def __init__(self, ff, units, density=25, seed=1, positive=True):
        """
        AuPlaneS core

        :param center_particle: type of the center particle
        :param units: the plane is made of units x units array of atoms
        :param positive: when true the sulfurs are built in the positive x direction from the Au plane
        :return: AuPlaneS object
        """

        # nanoparticle name
        if positive:
           plus = '+'
        else:
            plus = ''

        np.random.seed(seed)

        cells = [[0, 0, 0]]
        dist = 3.96
        y = np.multiply(range(-units,units+1,1), dist)
        z = np.multiply(range(-units,units+1,1), dist)

        for one in y:
            for two in z:
                cells.append([0,one,two])



        self.position = cp.deepcopy(cells)
        news = [[0, 1.98, 1.98], [0, -1.98, 1.98], [0, -1.98, -1.98], [0, 1.98, -1.98]]
        for cell in cells:
            for new in news:
                to_add = list(np.add(cell, new))
                if to_add not in self.position:
                    self.position.append(to_add)

        #centroid = np.average(self.position, axis=0)
        #self.position = np.subtract(self.position, centroid)
        total_max = np.max(self.position)
        to_remove = []
        for pos in self.position:
            if np.abs(np.max(pos) - total_max) < 10e-3:
                to_remove.append(pos)

        for tr in to_remove:
            self.position.remove(tr)
        golds = len(self.position)

        side = int(np.max([la.norm(c[1]) for c in self.position]) + 1.98)
        num = int((side*2)**2/density)
        #triangles = ((1+units)*2 -1)**2
        sulfurs = []
        chance = 3.92 * 4 /density
        length = 2 * (1+units) -1
        area = 3.96 ** 2 * (length + 1) ** 2
        nsulfurs = int(area/density)
        sulfurs = []
        spacing = np.sqrt(density)
        l = int(np.sqrt(nsulfurs))
        for i in range(1, l+1):
            for x in range(1,l+1):
                sulfurs.append([0, i *spacing, x*spacing])
        centroid = np.average(sulfurs, axis=0)
        sulfurs = np.subtract(sulfurs, centroid)

        if positive:
            sulfurs = [[2, s[1], s[2]] for s in sulfurs]
        else:
            sulfurs = [[-2, s[1], s[2]] for s in sulfurs]

        position = self.position

        super(AuPlaneS, self).__init__(ff, golds+1, 'AuPlane' + plus + str(units) + 'S')

        self.position = np.array([[0.0, 0.0, 0.0]] + list(position))

        self.mass[1:1+golds] = self.ff_reader.get_molecular_weight('Au')
        self.mass[1+golds:] = self.ff_reader.get_molecular_weight('S')
        self.mass[0] = np.sum(self.mass[1:])

        # types (Gold, Sulfur)
        self.types = ['_'+self.name, 'Au']
        self.typeid = ['_'+self.name] + ['Au']*golds
        self.graft_sites = sulfurs
        self.graft_num = len(sulfurs)

        self.moment_inertia[0] = np.diag(self.moment_of_inertia())

    def get_vector(self):
        """
        See Documentation in BasicSystemEntity
        """

        return np.cross(self.position[1] - self.position[0], self.position[2] - self.position[0])  # points perpendicular to the plane
