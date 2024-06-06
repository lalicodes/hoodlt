.. _AnalysisExample4:

This example illustrates how to obtain the tag numbers corresponding to a :mod:`FunctionalizedConfiguration`
A detailed explanation is given here :ref:`HOODLTConfigOrg`

Example 4: Conversion between FunctionalizedConfiguration and snapshot
======================================================================

.. code-block:: python

    import hoodlt.Data.Modelconfigurations.MapSnapshot as Ma
    from hoodlt.Analysis.Collect.ReInitHelper import reinit_config
    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds
    import hoodlt.Data.Modelconfigurations.Saver as sv

    name = 'cAu201S-Hydrocarbon-n11+cAu1072S-Hydrocarbon-n11_cMgzn2_l221_a48_ffDry-Ncs_uAngAmuEv_l_00p2'
    ff = 'dry-ncs'

    conf_nc = reinit_config(name)

    sim = SimulationWithBonds(sysfile=name, rcut=5, ff_name=ff, temp_in_kelvin=387)

    snap = sim.snap

    map_snap = Ma.MapSnapshot(conf_nc, snap)

    # get all the rigid center tags
    lst = map_snap.get_rigid_center_types()
    print(lst)

    # get tags of the different nanoparticles
    lst_ents = map_snap.list_of_entities
    # iterate over the entities
    for ents in lst_ents:
        txt = 'begins and ends at tag '
        print('name ', ents['tag-center-name'], txt, ents['tag-center'], ents['tag-end'])

    lst_tag = [176, 45816, 53300, 72613, 72614, 90469]
    typs = snap.particles.types
    print('list of all types', typs)
    for tag in lst_tag:
        print('tag ', tag, ' in object with core', map_snap.entity_from_tag(tag).core.get_name(), 'atom type ', typs[snap.particles.typeid[tag]])

    center_types = ['_cAu1072S', '_cAu201S']
    for cts in center_types:
        val=map_snap.rigid_body_params(cts)


This script provides a list (lst) in which each element in FunctionalizedObject
(particles, substrates or solvent) has a list entry. The elements of each list contain the
snapshot tags for each element. In this case, the system corresponds to 48 nanoparticles, so
therefore the list has 48 entries. Each entry is a dictionary, where for a particle

#. type: FunctionalizedParticle type: particles
#. center_name : name of the center (rigid body)
#. tag-center-name: tag of the center
#. tag-end: last tag for this particle
#. rigid: the tags that define rigid beads within this particles
#. flexible: tags associated with this particle that are flexible
