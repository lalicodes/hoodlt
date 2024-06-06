"""
:module: FccConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze an fcc unit cell configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ClusterConfigDataAbs import ClusterConfigDataAbs


class FccConfigData(ClusterConfigDataAbs):
    """
    Class which stores all the data for building and analyzing an fcc configuration of ncs
    """

    def __init__(self, dist_from_origin, nc_at_origin=True):
        """

        :param dist_from_origin: distance the cage ncs should be placed from the origin
        :param nc_at_origin: True if the position and bond list should include an nc at the origin
        """

        # define positions as scaled unit vectors
        positions = []
        t = 1/np.sqrt(2)
        unit_vectors = [(t, t, 0), (-t, t, 0), (-t, -t, 0), (t, -t, 0), (t, 0, t), (0, t, t), (-t, 0, t), (0, -t, t),
                        (t, 0, -t), (0, t, -t), (-t, 0, -t), (0, -t, -t)]

        if nc_at_origin:
            unit_vectors = [(0, 0, 0)] + unit_vectors

        for unit_vector in unit_vectors:
            positions.append([dist_from_origin*unit_vector[0],
                              dist_from_origin*unit_vector[1],
                              dist_from_origin*unit_vector[2]])

        nc_positions = np.array(positions)

        # define list of bonds
        num_cage_ncs = 12
        bond_list = []
        start_index = 0  # index to start when assigning cage bonds

        if nc_at_origin:
            start_index = 1
            for i in range(1, num_cage_ncs+1):
                bond_list.append((0, i, 0))

        bonds = [(0,4), (0,5), (0,8), (0,9),
                 (1,5), (1,6), (1,9), (1,10),
                 (2,6), (2,7), (2,10), (2,11),
                 (3,4), (3,7), (3,8), (3,11),
                 (5,4), (5,6), (7,4), (7,6), (9,8), (9,10), (11,8), (11,10)]

        for (i1, i2) in bonds:
            bond_list.append((i1+start_index, i2+start_index, start_index))

        # set ratios for analysis purposes
        bond_number_ratio = 2
        bond_length_ratio = 1

        super(FccConfigData, self).__init__(nc_positions, bond_list, bond_length_ratio, bond_number_ratio, num_cage_ncs)
