"""
:module: SimulationActions
:platform: Unix, Windows
:synopsis: Classes that operate with SimulationWithBonds

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2022
.. history:
"""


import numpy as np
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions


class ParticleDistance(object):

    def __init__(self, simulation, bond_params, local_snap):
        """
        constructor

        :param simulation: simulation object
        :param bond_params: bond parameters
        :param local_snap: hoomd snapshot
        """
        self._sim = simulation
        self._bond_params = bond_params
        self._local_snap = local_snap

    @property
    def all(self):
        """
        computes the distance for each and everyone of the bonds
        """
        dist = np.array([])
        snap = self._sim.state.get_snapshot()
        box = BoxFunctions(snap.configuration.box)
        for key in self._bond_params:
            ind_b = snap.bonds.types.index(key)
            indexes = [ind_x for ind_x, x in enumerate(snap.bonds.typeid) if x == ind_b]
            ind_mat = snap.bonds.group[indexes]
            pos1 = snap.particles.position[ind_mat[:, 0]]
            pos2 = snap.particles.position[ind_mat[:, 1]]
            dist = np.concatenate([dist, box.compute_distances(pos1, pos2)])
        return dist

    @property
    def average(self):
        """
        Computes the average distance for each bond type
        """
        dist = []
        snap = self._sim.state.get_snapshot()
        box = BoxFunctions(snap.configuration.box)
        for key in self._bond_params:
            ind_b = snap.bonds.types.index(key)
            indexes = [ind_x for ind_x, x in enumerate(snap.bonds.typeid) if x == ind_b]
            ind_mat = snap.bonds.group[indexes]
            pos1 = snap.particles.position[ind_mat[:, 0]]
            pos2 = snap.particles.position[ind_mat[:, 1]]
            dist.append(np.average(box.compute_distances(pos1, pos2)))
        return np.array(dist)

