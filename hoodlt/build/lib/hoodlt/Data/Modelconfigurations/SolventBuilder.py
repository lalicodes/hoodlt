"""
:module: SolventBuilder
:platform: Unix, Windows
:synopsis: Builds configurations of NCs

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> June 2022
.. history:
"""


import numpy as np
import numpy.linalg as la
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions


class SolventBuilder:
    """
    Class bookkeeping solvent locations
    """

    def __init__(self, box, grid_spacing, conf):
        """

        :param box: simulation box
        :param grid_spacing: spacing defining the grid
        :param conf: FunctionalizedConfiguration objects
        """

        particles = conf.particles
        substrates = conf.substrates
        solvent = conf.solvent

        self.actual_box = box
        self.spacing = grid_spacing

        self.box = BoxFunctions(self.actual_box)
        # maximum value
        self.max = 0.5*np.max(grid_spacing)
        # grid points
        self.points_grid = self.box.make_grid(self.spacing)
        # number of points in the grid
        self.points_number = self.points_grid.shape[0]
        # points in the grid that are available to place solvent
        self.points_avail = np.ones(self.points_number, dtype=bool)

        # find the available space for particles
        for part in particles:
            ct_pos = np.zeros([1, 3])
            ct_pos[0, :] = part.core.position[0]
            # compute distance of core surface to core center
            dist_core = self.box.compute_all_distances(part.core.position[1:, :], ct_pos)
            # compute max core radius
            # min_dist = np.amin(dist_core)
            max_dist = np.max(dist_core)
            dist = self.box.compute_all_distances(self.points_grid, ct_pos)
            self.points_avail[np.nonzero(dist < max_dist+self.max)[0]] = False
            for lig in part.ligands:
                dist = self.box.compute_all_distances(self.points_grid, lig.position)
                self.points_avail[np.nonzero(dist < self.max)[0]] = False
        # find the available space for substrates
        for subs in substrates:
            subs_vec = subs.get_vector()
            dist = np.abs(np.dot(self.points_grid - subs.position[0], subs_vec)) / la.norm(subs_vec)
            self.points_avail[np.nonzero(dist < self.max)] = False
        # find the available space for solvent
        for solv in solvent:
            dist = self.box.compute_all_distances(self.points_grid, solv.position)
            self.points_avail[np.nonzero(dist < self.max)[0]] = False

    def insert(self, solv):
        """
        adds one solvent molecule

        :param solv: solvent object
        """

        positions = solv.position
        max_dist = np.amax(la.norm(positions[np.newaxis, :, :]-positions[:, np.newaxis, :], axis=2))
        l_m = self.box.lsize()
        if max_dist > 0.5*np.min(l_m):
            txt = 'Size %1.2f solvent does not fit in box [%1.4f,%1.4f,%1.4f]' % (max_dist, l_m[0], l_m[1], l_m[2])
            raise ValueError(txt)
        solvent_placed = False
        # get all available indices
        ind_avail = np.nonzero(self.points_avail)
        for ind in ind_avail[0]:
            solv.shift(self.points_grid[ind])
            positions = solv.position
            dist = self.box.compute_all_distances(self.points_grid, positions)
            indices = np.unique(np.nonzero(dist < self.max)[0])
            if np.all(np.in1d(indices, ind_avail)):
                self.points_avail[indices] = False
                solvent_placed = True
                break
            else:
                solv.shift(-self.points_grid[ind])

        if not solvent_placed:
            raise ValueError('The solvent cannot be placed. Consider lowering the grid size')

    def add(self, solv, point):
        """
        adds one solvent molecule

        :param solv: solvent object
        :param point: point where to add the solvent
        """

        solv.shift(point)
        # get all available indices
        ind_avail = np.nonzero(self.points_number)
        positions = solv.position
        dist = self.box.compute_all_distances(self.points_grid, positions)
        indices = np.nonzero(dist < self.max)[0]
        if not set(indices).issubset(ind_avail):
            print('Warning: solvent maybe placed too tight')
        self.points_avail[indices] = False
