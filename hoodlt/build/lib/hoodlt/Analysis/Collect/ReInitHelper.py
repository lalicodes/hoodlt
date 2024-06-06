"""
:module: ReInitHelper
:platform: Unix, Windows
:synopsis: Functions which aid in reinitializing simulations through trajectory objects

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - units are now applied when copying information from the gsd, not after
..                  - Renamed the file
..                  - Added support for substrates
..                  - Added documentation
..
..                 Alex Travesset <trvsst@ameslab.gov> April 2022
..                  - made class compatible with the new unit definition
"""

import gsd.hoomd
import copy as cp
from hoodlt.Data.Modelconfigurations.Saver import load_config
from hoodlt.Data.Modelconfigurations.FunctionalizedParticle import FunctionalizedParticle
from hoodlt.Data.Modelsubstrates.SubstrateAbs import SubstrateAbs
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader


def reinit_config(pickle_name, gsd_name=None, with_restart=True):
    """
    Reinitializes a configuration from a (single frame) gsd file using the pickle file for some essential information.

    :param pickle_name: name of the pickle file, without the .pickle extension
    :param gsd_name: name of the gsd file, without the .gsd extension
    :param with_restart: include '_restart' in the name of the gsd file
    :return: FunctionalizedConfiguration object representing the first frame of the gsd file
    """

    # load the pickle and gsd file
    conf_old = load_config(pickle_name)
    units = ForceFieldReader(conf_old.ff_reader.name).get_units()
    if gsd_name is None:
        gsd_name = pickle_name
    sys_name = gsd_name
    if with_restart:
        sys_name += '_restart'
    sys = gsd.hoomd.open(sys_name + '.gsd', mode='r')

    # build the configuration
    conf_new = reinit_frame(conf_old, sys, 0, units)

    return conf_new


def reinit_frame(conf, sys, frame, units):
    """
    Reinitializes a given frame of a gsd file from the old configuration,

    :param conf: functionalizedConfiguration object. The result of unpickling a pickle file.
    :param sys: hoomd system object created from a gsd file
    :param frame: index of the frame in the system object
    :param units: units
    :return: a FunctionalizedConfiguration object which represents the state of the input system object at the 'frame'
    """

    if frame >= len(sys) or frame < 0:
        raise ValueError("frame out of bounds")

    snap = sys[int(frame)]
    conf_new = cp.deepcopy(conf)
    tag = 0

    for entity in conf_new.bondable_entities:
        if isinstance(entity, FunctionalizedParticle):
            tag = _copy_nc_chunk(entity, snap, tag, units)
        elif isinstance(entity, SubstrateAbs):  # at this point it must be a substrate
            n = entity.get_num_particles()
            tag = _copy_chunk(entity, 0, snap, tag, n, units)

    # solvents
    for solv in conf_new.solvent:
        n = solv.get_num_particles()
        tag = _copy_chunk(solv, 0, snap, tag, n, units)

    # box
    snap.configuration.box[:2] /= units.length_construction_to_simulation
    conf_new.create_box(snap.configuration.box)

    return conf_new


def _copy_nc_chunk(nc, snap, tag, units):
    """
    Copies a chunk of data from a snapshot starting at index tag to the nanocrystal being reinitialized

    :param nc: FunctionalizedParticle object
    :param snap: hoomd snapshot object
    :param tag: index in the snapshot object to begin copying
    :return: the tag after the copying is done
    """

    # core data
    core_count = nc.core.get_num_particles()
    tag = _copy_chunk(nc.core, 0, snap, tag, core_count, units)
    # ligand data, must copy graft sites twice
    for i, lig in enumerate(nc.ligands):
        n = lig.get_num_particles()
        tag = _copy_chunk(lig, 0, snap, tag, n, units)
        nc.core.graft_sites[i] = lig.position[0]  # update the graft positions on the core
    return tag


def _copy_grafter(core, core_tag, snap, snap_tag, units):
    """
    Copies the grafter in the snapshot object at index snap_tag to the core object at index core_tag

    :param core: the core object being reinitialized
    :param core_tag: the index in the array of core data to begin copying
    :param snap: the snapshot object that is being reinitialized
    :param snap_tag: the index in the snapshot object to begin copying
    :return: None
    """

    _copy_chunk(core, core_tag, snap, snap_tag, 1, units)


def _copy_chunk(basic_entity, basic_entity_tag, snap, snap_tag, size, units):
    """
    Copies a chunk of the snapshot object of length 'size' starting at index 'snap_tag' to the 'basic_entity' object at
    index 'basic_entity_tag'

    :param basic_entity: the core, ligand, or solvent object being reinitialized
    :param basic_entity_tag: index in the array of the basic entity's data to start copying
    :param snap: the snapshot object we are using for reinitialization
    :param snap_tag: the index in the snapshot object's data to start copying
    :param size: the length of the chunk of data to copy
    :return: the tag after copying is complete
    """

    quant_copy = ["position", "image", "velocity", "orientation"]
    quant_with_dimension = ["velocity", "acceleration", "charge", "mass", "moment_inertia"]  # position, diameter too
    for i in range(size):
        for quant in quant_copy:
            getattr(basic_entity, quant)[basic_entity_tag + i] = getattr(snap.particles, quant)[snap_tag+i]

    # remove the units from all the indexes you just copied
    basic_entity.position[basic_entity_tag: basic_entity_tag + size] /= units.length_construction_to_simulation
    basic_entity.diameter[basic_entity_tag: basic_entity_tag + size] /= units.length_construction_to_simulation
    txt = '_construction_to_simulation'
    for quant in quant_with_dimension:
        arr = getattr(basic_entity, quant)[basic_entity_tag: basic_entity_tag + size]
        arr /= getattr(units, quant+txt)
    tag = snap_tag + size

    return tag
