"""
:module: Squeeze
:platform: Unix, Windows
:synopsis: squeeze the NCs in the given configuration to the given radius and shape

.. moduleauthor:: Xun Zha <xzha@iastate.edu> Dec 2020
.. history:
..                Xun Zha <xzha@iastate.edu> June 2021
..                  - modified the class, now squeeze a list of FunctionalizedParticle objects,
..                    instead of a LatticeFunctionalizedConfiguration object
..
..                Alex Travesset <trvsst@ameslab.gov> April 2022
..                  - made the class to support hoomd v3
"""

import copy
import numpy as np
import numpy.linalg as la
import hoomd
from hoodlt.Data.Modelconfigurations.Saver import save_config
from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
from hoodlt.HOOMD.SimParameters import SimParameters
from hoodlt.Data.Forcefield.ScaledForceFieldReader import ScaledForceFieldReader
from hoodlt.Data.ProcessConfigurations.CreateWalls import CreateWalls


class Squeeze:
    def __init__(self, ff, list_nc, list_radius, wall_name='sphere', delta=None, dt=None, device='cpu'):
        """Takes a list of FunctionalizedParticle objects and squeezes each one into the desired radius and shape

        :param ff: forcefield name
        :param list_nc: list of NC(s) to squeeze
        :param list_radius: list of NC radius after squeezing
        :param wall_name: type of wall (sphere, cylinder, etc..)
        :param delta: fraction of the unit of distance by which the wall radius will be reduced at each time step
        :param dt: delta time step
        :param device: cpu or gpu
        :return:
        """

        # check parameters
        if not isinstance(list_nc, list):
            raise ValueError('Parameter "list_nc" needs to be a list of FunctionalizedParticle objects')
        if not isinstance(list_radius, list):
            raise ValueError('Parameter "list_radius" needs to be a list of NC radius')
        if len(list_nc) != len(list_radius):
            raise ValueError('Parameters "list_nc", "list_radius", "list_shape" should have the same length')

        self.ff = ff
        self.ff_reader = ScaledForceFieldReader(ff)
        self.units = self.ff_reader.get_units()

        self.list_nc = list_nc
        self.list_radius = list_radius

        # create walls
        self.wall_name = wall_name
        self.wall = CreateWalls(self.wall_name)
        # change how walls are build
        self.delta = delta
        if self.delta is None:
            self.delta = 0.1
        self.dt = dt
        if self.dt is None:
            self.dt = self.units.dt

        self.dev = device

    def squeeze(self):
        """Squeezes the NCs

        :return: list of squeezed NCs
        """
        new_list_nc = []
        for nc, radius in zip(self.list_nc, self.list_radius):
            new_nc = self._squeeze_single_nc(nc, radius)
            new_list_nc.append(new_nc)
        return new_list_nc

    def _squeeze_single_nc(self, nc, radius, sim_params=None, sigma_wall=None):
        """Squeezes a single NC

        :param nc: NC to squeeze
        :param radius: designated squeezed NC radius
        :Param sigma_wall: Lj value of the wall thickness
        :params sim_params: SimParam object
        :return:
        """

        # find the nc center and shift
        center_particle = nc.core.position[0]
        nc.shift(-center_particle)

        # determine sigma_wall
        if sigma_wall is None:
            sigma_wall = 3.0
        cut_wall = 1.1*sigma_wall

        # find initial NC radius and the increment radius to ensure that the final radius is reached exactly
        r_max = self._get_radius(nc)
        r_initial = r_max + sigma_wall
        r_new = r_initial
        r_final = radius
        num = (r_final-r_initial)/self.delta
        num_v = int(np.ceil(num))
        radius_increment = num*self.delta/num_v

        # build a single particle object with the nc info
        builder = ConfigurationBuilder()
        builder.add_reinit_nc(nc, np.zeros(3))
        builder.set_box(np.ones(3) * r_initial * 10)
        builder.set_alias('Squeeze' + '_' + self.wall_name)
        config = builder.get_configuration()
        # need to eliminate the _restart sufix from name
        name = save_config(config, device=self.dev)

        # HoomdClass parameters
        if sim_params is None:
            sim_params = SimParameters(500, 1000, [], 1000)

        simulation = HoomdSimulation(sysfile=name, rcut=5, ff_name=self.ff, write_data=False, translation=False)
        # adjust wall diameters
        while r_max > r_final:
            walls = self.wall.make_wall(r_new)
            lj_wall = hoomd.md.external.wall.LJ(walls=walls)
            lj_wall.params[nc.get_distinct_particle_types()] = {'sigma': sigma_wall, 'epsilon': 1, 'r_cut': cut_wall}
            simulation.potentials['wall'] = lj_wall
            simulation.run_hoomd(sim_params, dt=self.dt)
            r_max = self.wall.max_radius(simulation.sim.state.get_snapshot().particles.position)
            r_new = np.min([r_max+sigma_wall, r_new-radius_increment])

        # take snapshot of the final configuration
        snap = simulation.sim.state.get_snapshot()

        # express the nc in canonical coordinates, so that the core is in body frame
        return self._update_positions(nc, snap)

    @staticmethod
    def _get_radius(nc):
        """Gets the maximum radius of the given NC

        :param nc: FunctionalizedParticle object
        :return:
        """
        r0 = 0
        for lig in nc.ligands:
            dis = la.norm(lig.position, axis=1)
            r0 = np.max([r0, np.max(dis)])
        return r0

    @staticmethod
    def _update_positions(nc, snap):
        """Updates the positions of the given NC to the given positions

        :param nc: NC with positions to be updated
        :param snap: snapshot with new positions
        :return: NC with updated positions
        """

        # update only non-rigid particles
        index_non_rigid = np.argwhere(snap.particles.body < 0)
        # index of ligand atom in pos
        ind_l = nc.core.num_particles
        # update the positions of the ligand atoms
        for lig in nc.ligands:  # update the positions of the ligand atoms
            for i in range(lig.num_particles):
                if ind_l in index_non_rigid:
                    lig.position[i] = snap.particles.position[ind_l]
                ind_l += 1
        return nc
