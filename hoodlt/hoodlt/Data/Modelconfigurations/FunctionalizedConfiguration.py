"""
:module: FunctionalizedConfiguration
:platform: Unix, Windows
:synopsis: A representation of the state of a simulation using objects

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> April 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - removed an unecessary function
..                  - bonds can now be added and dumped between substrates and ncs
..                  - Moved methods that were in this class to Analysis.Analyze
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - calls SystemEntity constructor to set the forcefield reader
..                  - updated getting nc occupied volume to account for ncs with images
..                  - forcefield reader for the configuration is now created from the reader for the basic entities
..                  - altered remove and apply units slightly
..                  - made to inherit from SystemEntity
..                  - documentation
..                  - updated code for orienting/de-orienting the ncs in a configuration
..                  - removed some unnecessary methods
..                  - removed references to body_count variable when dumping
..                  - fixed issue with calculating masses on entities with rigid center particles
..                  - removed the '_' in front of many of the names of the methods
..                Tommy Waltmann <tomwalt@iastate.edu> May 2019
..                  - rewrote many of the methods to be simpler, easier to understand
..                  - added a few methods that get info needed to set up hoomd simulations
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                  - added special coulomb and lj functionality
..                Xun Zha <xzha@iastate.edu> July 2021
..                  - changed 'self.particles = nanos' to 'self.particles = copy.deepcopy(nanos)'
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                  - fixed a bug in ParticleToConfiguration
..                  - eliminated the dependence on hoomd classes
..                  - Made it consistent with HOOMD v3
"""
import copy
import copy as cp
import numpy as np
from hoodlt.Data.Modelconfigurations.SystemEntity import SystemEntity
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
from hoodlt.Data.Modelsubstrates.SubstrateAbs import SubstrateAbs


class FunctionalizedConfiguration(SystemEntity):
    """
    Defines a "snapshot" of a simulation. Methods that don't explicitly have documentation in this class have
    documentation in SystemEntity.
    """

    def __init__(self, nanos, positions, box_l, alias=None):
        """

        :param nanos: a list of FunctionalizedParticle objects
        :param positions: a list of positions corresponding to the locations of the center particles
        :param box_l: n-dimensional numpy array with optional len 3, 4, 6, defining the hoomd box parameters
        :param alias: A name given to the configuration by the user
        """
        if len(nanos) != len(positions):
            raise ValueError('must have the same number of particles as positions')
        for nc in nanos:
            if not np.allclose(nc.core.position[0], np.zeros(3)):
                raise ValueError('center cores must be placed at the origin')

        # initialize list of entities
        self.list_entities = ['particles', 'substrates', 'solvent']
        for entity in self.list_entities:
            setattr(self, entity, [])

        # keep distinct entities
        for entity in self.list_entities:
            setattr(self, 'distinct_'+entity, [])

        self.particles = cp.deepcopy(nanos)

        # keep all distinct functionalized particles
        list_names = []
        for ind, part in enumerate(self.particles):
            name = part.get_name()
            if name not in list_names:
                list_names.append(name)
                self.distinct_particles.append({name: copy.deepcopy(part)})

        # store list of ncs (at the shifted positions)
        self.bondable_entities = []
        for i in range(len(self.particles)):
            self.particles[i].shift(positions[i])
            self.bondable_entities.append(self.particles[i])

        # bond info for adding bonds between functionalized configurations
        self.bonds = []
        self.bonds_types = []
        self.bonds_dist = []

        # set the force field reader for the configuration
        if len(nanos) > 0:
            ff_reader = cp.deepcopy(nanos[0].core.ff_reader)  # this assumes all the ncs read from the same forcefield
        else:
            ff_reader = None

        super(FunctionalizedConfiguration, self).__init__(ff_reader)

        self.box = None
        self.create_box(box_l)

        self.alias = alias
        self.units_defined = 'construction'

    def create_box(self, box_l):
        """
        Creates the box from the input given from the user

        :param box_l: the box input from user
        :return: None
        """

        if isinstance(box_l, np.ndarray) or isinstance(box_l, list):
            if len(box_l) == 3:
                self.box = np.array([box_l[0], box_l[1], box_l[2], 0, 0, 0])
            elif len(box_l) == 4:
                self.box = np.array([box_l[0], box_l[0], box_l[0], box_l[1], box_l[2], box_l[3]])
            elif len(box_l) == 6:
                self.box = box_l
            else:
                print('box_l can only be 3-dim, 4-dim or 6-dim')

        elif isinstance(box_l, float) or isinstance(box_l, int):
            self.box = np.array([box_l, box_l, box_l, 0, 0, 0])

        else:
            raise ValueError('Box is not properly defined. Define simple cubic boxes as a float, and more complex boxes'
                             'with n-dimensional numpy arrays.')

    def get_num_particles(self):

        total = 0

        for entity_type in self.list_entities:
            for entity in getattr(self, entity_type):
                total += entity.get_num_particles()

        return total

    def add_to_snap(self, snap, tag):

        for entity in self.bondable_entities:
            tag = entity.add_to_snap(snap, tag)

        for solv in self.solvent:
            tag = solv.add_to_snap(snap, tag)

        # add data for bonds between nc centers
        self._add_ctr_bond_data_to_snap(snap)

        return tag

    def _add_ctr_bond_data_to_snap(self, snap):
        """
        Adds the information for the center-center bonds on this configuration to the snapshot

        :param snap: the snapshot object
        :return: None
        """

        # get all indices of rigid centers
        ind_rigid_nc_ctrs = self._get_indexes_of_bondable_centers()
        # get all bonds added so far, and include the new bonds that are built
        dist_b_types = self.get_distinct_bonds_data('bonds')
        # add new types if present
        for bnd in dist_b_types:
            if bnd not in snap.bonds.types:
                snap.bonds.types.append(bnd)
        # add the bonds
        for ind_b in range(len(self.bonds_types)):
            for bond in self.bonds[ind_b]:
                ind1, ind2 = bond
                snap.bonds.group.append((ind_rigid_nc_ctrs[ind1], ind_rigid_nc_ctrs[ind2]))
                snap.bonds.typeid.append(dist_b_types.index(self.bonds_types[ind_b]))

        snap.bonds.N = len(snap.bonds.group)

    def _get_indexes_of_bondable_centers(self):
        """
        Gets the indexes of the configuration rigid center particles on bondable entities within a gsd to be dumped.
        This is a helper method for the helper method that adds all the center-center bond data to a snapshot object.

        :return: a list of indexes in a gsd file to be dumped
        """

        result = np.zeros(len(self.bondable_entities), dtype=int)

        for i in range(1, len(result)):
            for j in range(i, len(result)):
                result[j] += self.bondable_entities[i-1].get_num_particles()

        return result

    def _get_box_for_dumping(self):

        return self.box

    def get_distinct_bonds_data(self, attr):
        """
        Helper method to get the distinct types of either particle, bond, angle, dihedral, or improper

        :param attr: the name of the type
        :return:
        """

        types = []

        for entity_type in self.list_entities:
            for entity in getattr(self, entity_type):
                for typ in getattr(entity, 'get_distinct_bonds_data')(attr):
                    if typ not in types:
                        types.append(typ)

        # need to also account for CTR-CTR bond types
        if attr == 'bonds':
            for typ in self.bonds_types:
                if typ not in types:
                    types.append(typ)

        return types

    def rotate(self, quat, indices=None):
        """

        :param quat: list of quaternions
        :param indices: indices of the particles to rotate
        :return: None
        """

        if indices is None:
            lst_quat = len(self.particles)*[quat[0]]
            indices = range(len(self.particles))
        else:
            if len(quat) != len(indices):
                raise ValueError('the number of quaternions is inconsistent with the number of particles to rotate')
            lst_quat = quat

        for ind, qt in enumerate(lst_quat):
            self.particles[indices[ind]].rotate(qt)

    def get_name(self):
        """
        provides the name for the object

        :return None:
        """

        name = ''

        # append core info, if core is present
        if len(self.particles) > 0:
            name += self._build_string(self.particles)
        elif len(self.substrates) > 0:
            name += 'substrate'
        else:
            name += 'pure_solvent'

        # append alias information, if it exists
        if self.alias is not None:
            name += ('_c' + self.alias)

        # append solvent info, if it exists
        if len(self.solvent) > 0:
            name += '_s' + self._build_string(self.solvent)

        # append substrate info
        if len(self.substrates) > 0:
            name += '_b' + self._build_string(self.substrates)

        # append forcefield info
        name += '_ff' + self.ff_reader.name.title()

        return name

    @staticmethod
    def _build_string(list_entities):
        """
        Builds a string from the distinct elements in the list by calling their get_name() method

        :param list: the list of SystemEntities to generate a string from
        :return: a string giving details about what is in the list
        """

        lst_all = [ent.get_name() for ent in list_entities]
        lst = []
        [lst.append(el) for el in lst_all if el not in lst]

        return '+'.join(lst)

    def apply_units(self):
        """
        Convert the configuration into simulation units
        """

        if self.units_defined == 'simulation':
            return

        units = ForceFieldReader(self.ff_reader.name).get_units()

        for entity_type in self.list_entities:
            for entity in getattr(self, entity_type):
                entity.apply_units()

        self.box = self.box*units.length_construction_to_simulation

        self.units_defined = 'simulation'

    def nullify(self):
        """
        Place all entities to their default values. This is used in pickling the object to ensure that the initial
        pickle has no data.
        """

        for entity_type in self.list_entities:
            lst_distinct = getattr(self, 'distinct_' + entity_type)
            lst = [list(val.keys())[0] for val in lst_distinct]
            lst_all = getattr(self, entity_type)
            for ind, entity in enumerate(lst_all):
                indx = lst.index(entity.get_name())
                getattr(self, entity_type)[ind] = cp.deepcopy(lst_distinct[indx][lst[indx]])
        # ensure that bondable entities point to the same object
        ind_p = 0
        ind_s = 0
        for ind, entity in enumerate(self.bondable_entities):
            if isinstance(entity, SubstrateAbs):
                obj = self.substrates[ind_s]
                ind_s += 1
            else:
                obj = self.particles[ind_p]
                ind_p += 1
            self.bondable_entities[ind] = obj

        self.units_defined = 'construction'


class ParticleToConfiguration(FunctionalizedConfiguration):
    """
    Converts a single particle to a FunctionalizedConfiguration object. Called by FunctionalizedParticle.Configuration()
    """

    def __init__(self, part):
        """
        This initializer should never actually be called; Use PairConfiguration or another class

        :param part: the FunctionalizedParticle object
        """

        # shift to the [0,0,0] center
        part.shift(-part.core.position[0])

        super(ParticleToConfiguration, self).__init__([part], [[0, 0, 0]], part._get_box_for_dumping())
