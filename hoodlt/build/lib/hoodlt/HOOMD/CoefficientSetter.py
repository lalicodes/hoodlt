"""
:module: CoefficientSetter
:platform: Unix, Windows
:synopsis: Defines helper functions that read directly from the force field

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, April 2019
                    - Rewrote every function to use the new ForcefieldReader/GsdParser classes
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - fixed typo in fslj and lj
..                  - added improper_harmonic function
..                  - nonbonded parameter functions now set a pairwise rcut instead of using system-wide rcut
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                  - added special coulomb and lj functionality
..                Jianshe Xia <xiajs6075@iccas.ac.cn> July 2021
..                  - added the table force field for bond, angle and dihedral
..                Alex Travesset <trvsst@ameslab.gov>, March 2022
..                  - Made it compatible with hoomd v3
"""
from hoodlt.Data.Modelconfigurations.Saver import load_config
from hoodlt.Data.Modelconfigurations.MapSnapshot import MapSnapshot


def set_rigid_bodies(rigid_object, snap, sysfile):
    """
    Sets the rigid body parameters for a hoomd simulation

    :param rigid_object: hoomd.md.Rigid object
    :param snap: snapshot
    :param sysfile: stystem file
    :return: None
    """

    conf = load_config(sysfile)
    map_snap_conf = MapSnapshot(conf, snap)
    list_center_types = map_snap_conf.get_rigid_center_types()

    for typ in list_center_types:
        rigid_object.body[typ] = map_snap_conf.rigid_body_params(typ)


def getattr_lowercase(st1, st2):
    """
    Converts lower case to Capitalized or Upper

    :param: string
    :param: string
    :return: an object
    """
    try:
        val = getattr(st1, st2.capitalize())
    except AttributeError:
        val = getattr(st1, st2.upper())

    return val


def set_electrostatic(list_electrostatic, rcut):
    """
    sets electrostatic parameters

    :param list_electrostatic: electrostatic parameter
    :param rcut: cut-off
    :return: dict
    """

    if list_electrostatic is None:
        list_electrostatic = [16, 4, rcut, 0]
        print('warning: default electrostatic cut-off is ', rcut, 'in units of simulation length')
    grid_linear, order_el, rcut_el, alpha_el = list_electrostatic
    grd = (grid_linear, grid_linear, grid_linear)

    dict_e = {'resolution': grd, 'order': order_el, 'r_cut': rcut_el, 'alpha': alpha_el}

    return dict_e


def set_nonbonded(hm_obj, nl, pair, snap, ff_reader, rcut):
    """
    sets non-bonded interactions

    :param hm_obj: hoomd pair object
    :param nl: neighborlist object
    :param pair: pair potential
    :param snap: hoomd snap
    :param ff_reader: force field reader
    :param rcut: cut-off
    """

    p_obj = getattr_lowercase(hm_obj, pair)(nlist=nl)
    l_part = snap.particles.types
    num_types = len(l_part)
    for ind1 in range(num_types):
        for ind2 in range(num_types):
            getattr(p_obj, 'params')[l_part[ind1], l_part[ind2]] = \
                ff_reader.get_non_bonded(pair, l_part[ind1], l_part[ind2])
            getattr(p_obj, 'r_cut')[l_part[ind1], l_part[ind2]] = \
                ff_reader.get_non_bonded_rcut(pair, l_part[ind1], l_part[ind2], rcut)

    return p_obj


def set_bonded(attr, hm_obj, state, potential, ff_reader):
    """
        sets bonded interactions

        :param attr: attribute, bonds, angles, etc...
        :param hm_obj: hoomd bonds, angle, etc.. object
        :param state: and snap.attr object
        :param potential: potential type: harmonic, etc
        :param ff_reader: force field reader
    """

    if potential == 'Table':
        width = int(ff_reader.get_all_width(attr))
        p_obj = getattr(hm_obj, potential)(width)
    else:
        p_obj = getattr_lowercase(hm_obj, potential)()  # atr = hoomd.md.bond.harmonic()
    for typ in getattr(state, 'types'):
        getattr(p_obj, 'params')[typ] = ff_reader.get_potentials_params(attr, potential, typ)

    return p_obj


def set_special(attr, hm_obj, state, potential, ff_reader, rcut):
    """
        sets bonded interactions

        :param attr: attribute, bonds, angles, etc...
        :param hm_obj: hoomd bonds, angle, etc.. object
        :param state: hoomd snap.state
        :param potential: potential type: harmonic, etc
        :param ff_reader: force field reader
        :param rcut: cutt off
    """

    p_obj = getattr_lowercase(hm_obj, potential)()
    for typ in getattr(state, 'types'):
        dict_v = ff_reader.get_potentials_params(attr, potential, typ)
        if 'defined_by_potential' in dict_v:
            pot_name = dict_v['defined_by_potential']
            alpha = dict_v['alpha']
            p1, p2 = _get_p_vals(typ)
            getattr(p_obj, 'params')[typ] = ff_reader.get_non_bonded(pot_name, p1, p2, alpha)
            getattr(p_obj, 'r_cut')[typ] = ff_reader.get_non_bonded_rcut(pot_name, p1, p2, rcut)
        else:
            getattr(p_obj, 'params')[typ] = dict_v
            getattr(p_obj, 'r_cut')[typ] = rcut

    return p_obj


def _get_p_vals(name):
    """
    gets the particle type for the first pair

    :param name: string
    :return : first particle
    """
    # we only need the first and last characters
    new_name = name.split('-')

    return new_name[0], new_name[-1]
