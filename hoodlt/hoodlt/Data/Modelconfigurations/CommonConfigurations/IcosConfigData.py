"""
:module: IcosConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze an icosahedral configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ClusterConfigDataAbs import ClusterConfigDataAbs


class IcosConfigData(ClusterConfigDataAbs):

    def __init__(self, dist_from_origin, nc_at_origin=True):
        """

        :param dist_from_origin: distance the cage ncs should be placed from the origin
        :param nc_at_origin: True if the position and bond list should include an nc at the origin
        """

        # define positions as scaled unit vectors
        positions = []

        # unit vectors for icos configuration, doing a converson from spherical coordinates
        unit_vectors = [(0, 0, 1), (0, 0, -1)]
        for i in range(10):
            theta = np.pi / 2 - (np.power(-1, i)) * np.arctan(1 / 2)
            phi = np.pi / 5 * i
            unit_vectors.append((np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)))

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

        bonds = [(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 2),
                 (2, 4), (4, 6), (6, 8), (8, 10), (10, 2),
                 (3, 5), (5, 7), (7, 9), (9, 11), (11, 3),
                 (0, 2), (0, 4), (0, 6), (0, 8), (0, 10),
                 (1, 3), (1, 5), (1, 7), (1, 9), (1, 11)]

        for (i1, i2) in bonds:
            bond_list.append((i1+start_index, i2+start_index, start_index))

        # set ratios for analysis purposes
        bond_number_ratio = 5/2
        bond_length_ratio = np.sqrt(2*(1-(1/np.sqrt(5))))

        super(IcosConfigData, self).__init__(nc_positions, bond_list, bond_length_ratio,
                                             bond_number_ratio, num_cage_ncs)
