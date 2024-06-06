"""
:module: CreateFilters
:platform: Unix, Windows
:synopsis: Defines filters to be used to calculate quantities

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> August 2022
.. history:
"""

import hoomd


def create_filter_single_entity(map_f_snap, include_entities='all'):
    """
    Creates a filter for each single entity

    :param map_f_snap: mapping between a functionalized configuration and a snapshot
    :param include_entities: include only the given entities
    :return: list of filters
    """

    list_ents = map_f_snap.list_of_entities
    list_filters = []

    if include_entities == 'all':
        list_single_ents = list(range(len(list_ents)))
    else:
        list_single_ents = include_entities

    for ind, ents in enumerate(list_ents):
        if ind in list_single_ents:
            tag_ini = ents['tag-center']
            # the end-tag is part of the structure
            tag_end = ents['tag-end']+1
            tags = list(range(tag_ini, tag_end))
            list_filters.append(hoomd.filter.Tags(tags))

    return list_filters
