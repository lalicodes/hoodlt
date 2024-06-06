"""
:module: BasicSystemEntity
:platform: Unix, Windows
:synopsis: Abstract class that can represent any basic system entity (ligand, core, solvent)
.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - removed unecessary functions
..                  - calls SystemEntity constructor to set the forcefield reader
..                  - altered remove and apply units slightly
..                  - made to inherit from SystemEntity
..                  - quick fix for calculating correct mass and center of mass
..                  - added impropers
..                  - added documentation
..                  - added align to this class
..                  - updated code for de-orienting/orienting system entities
..                  - removed some unnecessary methods
..                  - removed references to body_count variable when dumping
..                  - fixed issue with calculating masses on entities with rigid center particles
..                  - removed the '_' in front of many of the names of the methods
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                   - added special coulomb and lj functionality
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                   - Made the class compatible with hoomd v3.0
..                   - Eliminated almost all references to explicit snap parameters
..                   - Made use of class decorators
..                Alex Travesset <trvsst@ameslab.gov> July 2022
..                   - fixed an important bug in regards to adding objects that include different types
"""


import numpy as np
import numpy.linalg as la
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Data.Modelconfigurations.SystemEntity import SystemEntity
from hoodlt.Data.Modelconfigurations.SnapshotH import SnapshotDefault
from hoodlt.Groups.Quat import Quat


def get_attrs(attr):
    """
    Returns the name of the function that returns all attrs

    :param attr: attribute (bond, angle, etc..)
    """

    return 'get_' + attr


def get_attr_types(attr):
    """
    Returns the name of the function that returns all attribute types

    :param attr: attribute (bond, angle, etc..)
    """

    return 'get_' + attr[:-1] + '_types'


def get_distinct_attrs():
    """
    Returns the name of the function that returns all distinct attribute type
    """

    return 'get_distinct_bonds_data'


def decor_rm_duplicates():
    """
    function decorator that removes duplicates from a list
    """
    def _decor(self, attr):
        lst = getattr(self, get_attr_types(attr))()
        lst_dist = []
        for val in lst:
            if val not in lst_dist:
                lst_dist.append(val)
        return lst_dist
    return _decor


def add_to_snap_bond(snap, tag, attr):
    """
    function adding bond attributes to a snap

    :param snap: HOOMD snapshot object
    :param tag: tag of the given particle
    :param attr: attribute (bond, angle, etc..
    """
    def _add(self):
        attr_types = getattr(self, get_attr_types(attr))()
        dist_attr_types = getattr(self, get_distinct_attrs())(attr)
        attribs = getattr(self, get_attrs(attr))()
        snap_attr = getattr(snap, attr)
        typs = getattr(snap_attr, 'types')
        for n_typ in dist_attr_types:
            if n_typ not in typs:
                typs.append(n_typ)
        for ind in range(len(attribs)):
            getattr(snap_attr, 'group').append(list(np.add(attribs[ind], tag)))
            getattr(snap_attr, 'typeid').append(typs.index(attr_types[ind]))
        setattr(snap_attr, 'N', len(getattr(snap_attr, 'group')))
    return _add


def add_to_snap_constraint(snap, tag, attr):
    """
    function adding constraint attributes to a snap

    :param snap: HOOMD snapshot object
    :param tag: tag of the given particle
    :param attr: attribute (bond, angle, etc..
    """
    def _add(self):
        attr_values = getattr(self, get_attr_types(attr))()
        attribs = getattr(self, get_attrs(attr))()
        snap_attr = getattr(snap, attr)
        for ind in range(len(attribs)):
            getattr(snap_attr, 'group').append(list(np.add(attribs[ind], tag)))
        
        setattr(snap_attr, 'value', getattr(snap_attr, 'value')+attr_values)
        setattr(snap_attr, 'N', len(getattr(snap_attr, 'group')))
    return _add


def define_bond_constraint_functions(cls):
    """
    class decorator that introduces bond and constraint functions
    """

    snap_d = SnapshotDefault()
    lst = snap_d.bonds_data + snap_d.constraints_data

    for attr in lst:
        setattr(cls, get_attrs(attr), lambda x: [])
        setattr(cls, get_attr_types(attr), lambda x: [])
    setattr(cls, get_distinct_attrs(), decor_rm_duplicates())

    return cls


@define_bond_constraint_functions
class BasicSystemEntity(SystemEntity):
    """
    Abstract class which can represent any "basic" system entity (i.e. any system entity which is not broken down
    further by storing instances of other predefined classes). Methods without documentation in this class have
    documentation in SystemEntity.
    """

    def __init__(self, forcefield, num_particles, name):
        """

        :param forcefield: the name of the forcefield to be used to build this entity
        :param num_particles: the number of particles that make up this entity
        :param name: the name of the entity
        """
        ff = ForceFieldReader(forcefield)

        super(BasicSystemEntity, self).__init__(ff)

        self.name = name
        self.num_particles = num_particles

        # set all quantities to be stored to their defaults, subclasses will override some of these quantities
        self.quantities = SnapshotDefault().particles_default

        for quant, default in self.quantities.items():
            setattr(self, quant, np.array([default]*num_particles))

        # inheriting classes must override this
        self.typeid = []

    def get_mass(self):
        """
        Gets the true mass of this entity. The mass of the rigid center particles is equal to the sum of the masses of
        its constituent particles, which creates errors in getting the mass when it is obtained by np.sum(self.mass).
        To get the correct value for the mass of the entity, use this method

        :return: the true mass of the entity
        """

        mass = 0
        for ind in self._get_indexes_rigid_centers_and_nonrigid_particles():
            mass += self.mass[ind]

        return mass

    def _get_indexes_rigid_centers_and_nonrigid_particles(self):
        """
        Gets the indexes of all the particles which are either rigid centers or nonrigid particles

        :return: a list of indexes of the particles which are not rigid centers
        """

        indexes = []
        for i in range(self.num_particles):
            if self.typeid[i][0] == 'c' or self.typeid[i][0] == '_':
                indexes.append(i)
            elif self.body[i] == -1:
                indexes.append(i)

        return indexes

    def _get_indexes_rigid_center_particles(self):
        """
        Gets the indexes of all of the particles which are rigid center particles

        :return: a list of indexes of particles which are rigid center particles
        """

        indexes = []
        for i in range(self.num_particles):
            if self.typeid[i][0] == 'c' or self.typeid[i][0] == '_':
                indexes.append(i)

        return indexes

    def get_distinct_particle_types(self):

        return self.types

    def get_name(self):
        """
        Gets the name of this entity
        :return: the name of this entity, a string
        """

        return self.name

    def get_num_particles(self):
        """
        Gets the number of particles on this entity

        :return: an integer
        """

        return self.num_particles

    def add_to_snap(self, snap, tag):
        d_snap = SnapshotDefault()
        for bnd in d_snap.bonds_data:
            add_to_snap_bond(snap, tag, bnd)(self)
        for cnd in d_snap.constraints_data:
            add_to_snap_constraint(snap, tag, cnd)(self)

        tag = self.add_particle_data_to_snap_range(self.num_particles, snap, tag)

        return tag

    def add_particle_data_to_snap_range(self, size, snap, tag):
        """
        Adds the first 'size' particles of this entity's particle data to the snapshot

        :param size: number of particles out of this entity to add to the snapshot
        :param snap: the snapshot object
        :param tag: the index in the snapshot to begin adding the entity's data at
        :return: the index in the snapshot immediately after the new data for this entity has been added
        """

        tag_start = tag
        get_dist_p_types = self.get_distinct_particle_types()
        # check if the types are already present
        for typ in get_dist_p_types:
            if typ in snap.particles.types:
                continue
            snap.particles.types.append(typ)

        dist_p_types = snap.particles.types
        for i in range(size):
            for quant in self.quantities:
                if quant != 'typeid':
                    val = getattr(self, quant)[i]
                    getattr(snap.particles, quant)[tag] = val
                else:
                    snap.particles.typeid[tag] = dist_p_types.index(self.typeid[i])

            # add body data to snap
            if self.body[i] > -1:
                snap.particles.body[tag] = tag_start + self.body[i]

            tag += 1

        return tag

    def center_of_mass(self):
        """
        Calculates the center of mass of this entity

        :return: a numpy array containing the position of the center of mass
        """
        r_cm = np.zeros(3)
        for i in self._get_indexes_rigid_centers_and_nonrigid_particles():
            for j in range(3):
                r_cm[j] += self.mass[i]*self.position[i, j]

        return r_cm / self.get_mass()

    def shift(self, vector):
        """
        Shifts the position of each particle on this entity by a given vector

        :param vector: translation vector
        :return: void
        """

        self.position[:] = self.position[:] + vector

    def rotate(self, relative_center, quat):
        """
        Rotates the basic entity according to the center and the rotation defined by the quaternion

        :param relative_center: position vector to use as the origin for the rotation
        :param quat: quaternion (as a tuple)
        """

        qt = Quat()
        for ind in range(len(self.position)):
            self.position[ind] = relative_center + qt.rotation(quat, self.position[ind]-relative_center)

    def align(self, relative_center, vector):
        """
        Aligns the vector obtained by calling the get_vector() method with the input vector through a rotation which
        uses relative center as its origin

        :return: void
        """

        self.align_using_vectors(relative_center, self.get_vector(), vector)

    def align_using_vectors(self, relative_center, original_vector, vector_to_rotate_to):
        """
        Rotates the entity in the exact same way that original_vector must be rotated so that it points in the same
        direction as vector_to_rotate_to

        :param relative_center: position to be used as the origin of the rotation
        :param original_vector: vector that will be used as the start reference for the rotation of the entity
        :param vector_to_rotate_to: vector that will serve as end reference for the rotation of the entity
        :return: None
        """

        q = Quat().vector_vector_rotation(original_vector, vector_to_rotate_to)

        self.rotate(relative_center, q)

    def get_vector(self):
        """
        Returns the reference vector used for orienting, de-orienting, and aligning this entity. It should be defined
        in terms of some positions on the entity, do not hard code this vector.

        :return: A vector (array length 3) which is used to do geometric rotations on this entity
        """

        raise ValueError("Must Override get_vector() method of BasicSystemEntity")

    def __eq__(self, other):
        """
        Determines whether this entity is the same as another entity

        :param other: the other BasicSystemEntity
        :return: True if this entity is the same as another entity, False otherwise
        """

        v = self.name == other.name and self.num_particles == other.num_particles \
            and self.ff_reader.name == other.ff_reader.name
        return v

    def scale(self, val):
        """
        takes all the positions and multiples them by the val; essentially rescaling the system

        :param val: scaling factor
        :return: void
        """

        self.position = self.position * val

    def moment_of_inertia(self, ctol=1e-3):
        """
        returns the moment of inertia. Intended to only be used on entities which are a single rigid body, with a rigid
        center particle located at the center of mass

        :param: ctol tolerance for center of mass
        :return: moment of inertia matrix (3,3)
        """

        # Make sure center of mass is at origin
        if la.norm(self.center_of_mass()) > ctol:
            y = self.center_of_mass()
            self.shift(-1*y)
            x = self.moment_of_inertia()
            self.shift(y)
            return x

        # Get matrix
        masses = [[mass] for mass in self.mass]
        i_diag = np.sum(masses*self.position**2)

        loc_inertia = i_diag*np.eye(3, 3) - np.tensordot(masses*self.position, self.position, axes=[0, 0])

        return loc_inertia

    def apply_units(self):

        units = ForceFieldReader(self.ff_reader.name).get_units()
        d_snap = SnapshotDefault()

        for key, value in d_snap.particles_units.items():
            if value is None:
                continue
            unit = getattr(units, value+'_construction_to_simulation')
            quant = getattr(self, key)
            setattr(self, key, quant*unit)
