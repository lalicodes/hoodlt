"""
:module: Saver
:platform: Unix, Windows
:synopsis: Saves FunctionalizedConfiguration objects

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> November 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - now the saver writes to a gsd only after converting to dimensionless units
..                  - Updated documentation
..
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                  - Modified the class to make it compatible with hoomd v3
"""

import pickle
import json


def save_config(f_object, name=None, add_restart=True, device='cpu'):
    """
    Writes a FunctionalizedConfiguration object to a gsd file and pickle file after converting to simulation units

    :param f_object: FunctionalizedConfiguration object
    :param name: an optional name to use for saving the configuration
    :param add_restart: end the file with the string _restart
    :device: gpu or cpu
    :return: the name of the gsd/pickle file the configuration was written to, without the .gsd extension
    """

    # set the save name
    if name is not None:
        save_name = name
    else:
        save_name = f_object.get_name()

    save_name_with_units = save_name + '_u' + f_object.ff_reader.get_units().name

    # dump information for bonds
    b_types = f_object.bonds_types
    b_dist = f_object.bonds_dist
    b_lamda = [1.0] * len(b_dist)
    dict_data = {'name': save_name_with_units, 'types': b_types, 'dist': b_dist, 'lamda': b_lamda}

    save_name_dist = save_name_with_units + '_bonds.json'
    with open(save_name_dist, 'w') as fp:
        json.dump(dict_data, fp)

    save_name_restart = save_name_with_units
    if add_restart:
        save_name_restart += '_restart'

    # apply units
    f_object.apply_units()

    # write to gsd
    f_object.write_gsd(save_name_restart, dev=device)

    with open(save_name_with_units + '.pickle', 'wb') as pickle_out:
        # erase all data in the object, this way pickle does not contain state information
        f_object.nullify()
        pickle.dump(f_object, pickle_out)

    return save_name_with_units


def load_config(name):
    """
    Unpickles a configuration that was saved using save_config()

    :param name: name of the pickle file, without the .pickle extension
    :return: a FunctionalizedConfiguration object
    """

    name_pickle = name + '.pickle'

    with open(name_pickle, 'rb') as filename:
        f_object = pickle.load(filename)

    return f_object
