"""
:module: TetraConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze a tetrahedral configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ClusterConfigDataAbs import ClusterConfigDataAbs


class TetraConfigData(ClusterConfigDataAbs):
    """
    Class which stores all the data necessary to build and analyze a tetrahedral configuration of ncs
    """

    def __init__(self, dist_from_origin, nc_at_origin=True):
        """

        :param dist_from_origin: distance the cage ncs should be placed from the origin
        :param nc_at_origin: True if the position and bond list should include an nc at the origin
        """

        # define positions as scaled unit vectors
        positions = []
        unit_vectors = [( np.sqrt(2/3), 0,-np.sqrt(1/3)),
                        (-np.sqrt(2/3), 0,-np.sqrt(1/3)),
                        ( 0, np.sqrt(2/3), np.sqrt(1/3)),
                        ( 0,-np.sqrt(2/3), np.sqrt(1/3))]

        if nc_at_origin:
            unit_vectors = [(0, 0, 0)] + unit_vectors

        for unit_vector in unit_vectors:
            positions.append([dist_from_origin*unit_vector[0],
                              dist_from_origin*unit_vector[1],
                              dist_from_origin*unit_vector[2]])

        nc_positions = np.array(positions)

        # define list of bonds
        num_cage_ncs = 4
        bond_list = []
        start_index = 0  # index to start when assigning cage bonds
        
        if nc_at_origin:
            start_index = 1
            for i in range(1, num_cage_ncs+1):
                bond_list.append((0, i, 0))

        for i in range(start_index, start_index+(num_cage_ncs-1)):
            for j in range(i+1, start_index+num_cage_ncs):
                bond_list.append((i, j, start_index))

        # set ratios for analysis purposes
        bond_number_ratio = 3/2
        bond_length_ratio = 2 * np.sqrt(2/3)

        super(TetraConfigData, self).__init__(nc_positions, bond_list, bond_length_ratio, bond_number_ratio, num_cage_ncs)
