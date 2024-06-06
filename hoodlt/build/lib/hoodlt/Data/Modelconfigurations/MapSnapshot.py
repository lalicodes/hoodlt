"""
:module: MapSnapshot
:platform: Unix, Windows
:synopsis: Functions which aid in reinitializing simulations through trajectory objects

.. moduleauthor:: Alex Travesset <trvsst@amelab.gov> April 2022
.. history:
"""
import numpy as np

from hoodlt.Data.Modelconfigurations.FunctionalizedConfiguration import FunctionalizedConfiguration
from hoodlt.Data.Modelconfigurations.FunctionalizedParticle import FunctionalizedParticle
from hoodlt.Data.Modelsubstrates.SubstrateAbs import SubstrateAbs
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class MapSnapshot:
    """
    Maps FunctionalizedConfiguration to snapshot and viceversA
    """

    def __init__(self, conf, snap):
        """
        Constructor

        :param conf: functionalizedConfiguration object. The result of unpickling a pickle file.
        :param snap: hoomd snapshot
        :return: list of tags for each entity
        """

        self.conf = conf
        self.snap = snap

        entities = FunctionalizedConfiguration([], [], [-1, -1, -1]).list_entities
        self.entities = entities

        self.list_of_entities = []
        snap_current = 0
        parts = snap.particles
        # FunctionalizedParticles and substrates
        for ind, entity in enumerate(self.conf.bondable_entities):
            pt = {'type': entities[0], 'tag-center': snap_current, 'position': ind}
            lst_rigid = [ind for ind in range(parts.N) if (parts.body[ind] == snap_current and ind != snap_current)]
            pt['rigid'] = lst_rigid
            snap_current = snap_current + entity.get_num_particles()
            pt['tag-end'] = snap_current - 1
            all_list = range(pt['tag-center'], snap_current)
            pt['flexible'] = list(set(all_list) - (set(pt['rigid']) | {pt['tag-center']}))
            if isinstance(entity, FunctionalizedParticle):
                if len(lst_rigid) == 0:
                    raise ValueError('FunctionalizedParticle has no rigid bodies')
                pt['tag-center-name'] = entity.core.types[0]
            elif isinstance(entity, SubstrateAbs):  # at this point it is a substrate
                pt['tag-center-name'] = entity.types[0]
            else:
                raise ValueError('incorrect value in bondable entities')
            self.list_of_entities.append(pt)
        # solvent
        for ind, solv in enumerate(self.conf.solvent):
            sl = {'type': entities[2], 'tag-center': snap_current, 'position': ind}
            lst_rigid = [ind for ind in range(parts.N) if (parts.body[ind] == snap_current and ind != snap_current)]
            snap_current = snap_current + solv.get_num_particles()
            sl['tag-end'] = snap_current - 1
            if len(lst_rigid) == 0:
                sl['rigid'] = None
                sl['tag-center-name'] = None
                sl['flexible'] = list(range(sl['tag-center'], sl['tag-end']))
            else:
                sl['rigid'] = lst_rigid
                sl['tag-center-name'] = solv.types[0]
                sl['flexible'] = None
            self.list_of_entities.append(sl)

    def entity_from_tag(self, tag_number):
        """snap_current
        Returns the entity that contains the given tag_number

        :param tag_number: hoomd snapshot tag
        :return: entity (FunctionalizedParticle, Substrate, Solvent)
        """

        for lst in self.list_of_entities:
            if lst['tag-center'] <= tag_number <= lst['tag-end']:
                ind = lst['position']
                if lst['type'] == 'solvent':
                    return self.conf.solvent[ind]
                else:
                    return self.conf.bondable_entities[ind]

    def get_rigid_center_types(self):
        """
        Gets a list of distinct types of rigid centers in the system.

        :return: a list of distinct types of rigid centers in the system.
        """

        list_ctrs = [st for st in self.snap.particles.types if st.startswith('_')]

        return list(set(list_ctrs))

    def rigid_body_params(self, center_type):
        """
        Returns the parameters for the given rigid_body

        :param center_type: the name of the center type
        """

        dc = {'constituent_types': None, 'positions': None, 'orientations': None}

        # find the rigid body in the configuration
        list_entity = ['particles', 'substrates', 'solvent']
        for entity in list_entity:
            for dict_p in getattr(self.conf, 'distinct_'+entity):
                for key, part in dict_p.items():
                    if isinstance(part, FunctionalizedParticle):
                        if center_type == part.core.typeid[0]:
                            # center is in functionalized particle object
                            # check that the snapshot does contain that center
                            if center_type not in self.snap.particles.types:
                                raise ValueError('FunctionalizedConfiguration and snap are inconsistent')
                            # number of constituent particles
                            num_p = part.core.num_particles + len(part.core.graft_sites) - 1
                            # get the rigid body types
                            for dicts in self.list_of_entities:
                                if dicts['tag-center-name'] == center_type:
                                    indx_rigid = dicts['rigid']
                                    typs = self.snap.particles.types
                                    list_types = [typs[self.snap.particles.typeid[ind]] for ind in indx_rigid]
                                    list_charges = [self.snap.particles.charge[ind] for ind in indx_rigid]
                                    break
                            # check that the number of particles is consistent
                            if num_p != len(list_types):
                                raise ValueError('FunctionalizedConfiguration and snap are inconsistent in particle')
                            # get the rigid body positions in body frame
                            list_pos = np.zeros([num_p, 3])
                            list_pos[:part.core.num_particles-1, :] = part.core.position[1:, :]
                            list_pos[part.core.num_particles-1:, :] = part.core.graft_sites[:, :]
                            break
                    elif isinstance(part, SubstrateAbs):
                        if center_type == part.typeid[0]:
                            if center_type not in self.snap.particles.types:
                                raise ValueError('substrate and snap are inconsistent')
                            num_p = part.num_particles - 1
                            # get the rigid body types
                            for dicts in self.list_of_entities:
                                if dicts['tag-center-name'] == center_type:
                                    indx_rigid = dicts['rigid']
                                    typs = self.snap.particles.types
                                    list_types = [typs[self.snap.particles.typeid[ind]] for ind in indx_rigid]
                                    list_charges = [self.snap.particles.charge[ind] for ind in indx_rigid]
                                    break
                            if num_p != len(list_types):
                                raise ValueError('FunctionalizedConfiguration and snap are inconsistent in substrate')
                                # get the rigid body positions in body frame
                            list_pos = np.zeros([num_p, 3])
                            list_pos[:] = part.position[1:]
                            break
                    elif isinstance(part, SolventAbs):
                        if center_type == part.typeid[0]:
                            num_p = part.num_particles - 1
                            for dictv in self.list_of_entities:
                                if dictv['tag-center-name'] == center_type:
                                    indx_rigid = dictv['rigid']
                                    typs = self.snap.particles.types
                                    list_types = [typs[self.snap.particles.typeid[ind]] for ind in indx_rigid]
                                    list_charges = [self.snap.particles.charge[ind] for ind in indx_rigid]
                                    break
                            list_pos = np.zeros([num_p, 3])
                            list_pos[:num_p, :] = part.position[1:, :]
                            break
                    else:
                        raise ValueError('center provided is not a rigid body')

        dc['constituent_types'] = list_types
        dc['positions'] = list_pos
        dc['orientations'] = num_p*[(1, 0, 0, 0)]

        return dc
