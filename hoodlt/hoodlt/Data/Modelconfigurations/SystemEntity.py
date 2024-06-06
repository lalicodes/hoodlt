"""
:module: SystemEntity
:platform: Unix, Windows
:synopsis: Abstract class which has methods that each system entity (configuration, nanoparticle, core, ligand, solvent)
must implement

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - removed an unecessary function
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - now sets a forcefield reader in the constructor, since every entity is built from a forcefield
..                  - altered remove and apply units slightly
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                  - Added special coulomb and lj functionality
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                  - Made the class consistent with HOOMD v3 and beyond
"""

import numpy as np
import gsd.hoomd
import hoomd
from hoodlt.Data.Modelconfigurations.SnapshotH import SnapshotDefault
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions


class SystemEntity:
    """
    Abstract class made to handle dumping of every type of entity through abstract methods, which will be overridden by
    inheriting classes
    """
    def __init__(self, ff_reader):
        """

        :param ff_reader: Forcefield reader object
        """

        self.ff_reader = ff_reader

    def get_name(self):
        """
        Returns the name of this entity, used as the default name for dumping to gsd and xyz files.

        :return: String, which is the name of the entity
        """

        self._raise_override_error('get_name')

    def get_num_particles(self):
        """
        Gets the number of particles in this entity

        :return: the number of particles in this entity
        """

        self._raise_override_error('get_num_particles')

    def get_distinct_particles(self):
        """
        Gets all the distinct particle types in this entity

        :return: a list of unique particle types in the entity
        """

        self._raise_override_error('get_distinct_particles')

    def get_distinct_bonds_data(self, b_types):
        """
        Gets all the distinct bond types in this entity
        :param b_types: bond types, like bonds, angles, etc..

        :return: a list of unique bond types in the entity
        """
        self._raise_override_error('get_distinct_bonds_data ' + b_types)

    def get_distinct_constraints_data(self, c_types):
        """
        Gets all the distinct constraint types in this entity

        :param c_types: constraints
        :return: a list of the unique special coulomb types in the entity
        """

        self._raise_override_error('get_distinct_constraints_data ' + c_types)

    def add_to_snap(self, snap, tag):
        """
        Adds this configuration to the snapshot object starting at index 'tag'

        :param snap: the snapshot object
        :param tag: the index to begin adding the configuration's data at
        :return: the index in the snapshot immediately after all the new data for this configuration has been added
        """

        self._raise_override_error('add_to_snap')

    def apply_units(self):
        """
        Scales all quantities with dimension on this entity. All quantities with dimension go from simulation units
        to dimensionless units

        :return: None
        """

        self._raise_override_error('apply_units')

    def make_snapshot(self):
        """
        Creates a snapshot object for this entity

        :return: a hoomd snapshot object
        """

        snap = gsd.hoomd.Frame()

        num_particles = int(self.get_num_particles())
        d_snap = SnapshotDefault()

        setattr(snap.particles, d_snap.N, num_particles)
        setattr(snap.particles, d_snap.types, [])

        for quant, default in d_snap.particles_default.items():
            setattr(snap.particles, quant, np.array([default] * num_particles))

        for bnd in d_snap.bonds_data:
            dict_b = getattr(d_snap, bnd+'_default')
            for key, value in dict_b.items():
                setattr(getattr(snap, bnd), key, value)

        for cnd in d_snap.constraints_data:
            dict_c = getattr(d_snap, cnd + '_default')
            for key, value in dict_c.items():
                setattr(getattr(snap, cnd), key, value)

        setattr(snap.configuration, d_snap.box, self._get_box_for_dumping())

        # add all particle, bond, angle, dihedral, improper, and special pair info to the snapshot
        self.add_to_snap(snap, 0)

        # return snapshot
        return snap

    def write_gsd(self, filename=None, dev='cpu', sd=36):
        """
        writes the entity to a gsd file

        :param filename: the name of the gsd file to dump to, without the .gsd extension. Defaults to what is returned
        by the get_name() method
        :param dev: device to use
        :param sd: seed
        :return: snapshot
        """

        snap = self.make_snapshot()

        # enforce periodic boundary conditions
        box = BoxFunctions(snap.configuration.box)
        pos_periodic = box.wrap(snap.particles.position)
        pos_image = box.get_images(snap.particles.position)
        snap.particles.position = pos_periodic
        snap.particles.image = pos_image

        file_name = self.get_name() + '.gsd'
        if filename is not None:

            file_name = filename + '.gsd'

        if dev == 'gpu':
            dv = hoomd.device.GPU()
        else:
            dv = hoomd.device.CPU()

        sim = hoomd.Simulation(device=dv, seed=sd)
        sim.create_state_from_snapshot(snap)
        hoomd.write.GSD.write(state=sim.state, filename=file_name, mode='wb')
        
        return sim.state

    def _get_box_for_dumping(self):
        """
        Returns and ad-hoc hoomd box object for dumping this entity to a gsd file

        :return: a hoomd box object
        """

        self._raise_override_error('_get_box_for_dumping')

    def _raise_override_error(self, method_name):
        """
        Raises the value error if a given method is not overridden a an inheriting class.

        :param method_name: the name of the method that was not overridden by an inheriting class
        :return: Nonesnap.particles.N = int(self.get_num_particles())
        """
        msg = "Must override abstract method " + method_name + "() of SystemEntity"
        raise ValueError(msg)
