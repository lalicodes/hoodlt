"""
:module: PlanarConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze a planar configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ClusterConfigDataAbs import ClusterConfigDataAbs


class PlanarConfigData(ClusterConfigDataAbs):
    """
    Class which stores all the data necessary to build and analyze symmetric, planar configurations of ncs
    """

    def __init__(self, num_cage_ncs, dist_from_origin, nc_at_origin=True):
        """

        :param num_cage_ncs: the number of ncs in the cage
        :param dist_from_origin: distance the cage ncs should be placed from the origin
        :param nc_at_origin: True if the position and bond list should include an nc at the origin
        """

        if num_cage_ncs < 2:
            raise ValueError("num_cage_ncs must be greater than or equal to 2")

        # define positions as scaled unit vectors
        positions = []
        if nc_at_origin:
            positions.append((0, 0, 0))
        for j in range(num_cage_ncs):
            positions.append([dist_from_origin*np.cos(2*np.pi*j/num_cage_ncs),
                              dist_from_origin*np.sin(2*np.pi*j/num_cage_ncs), 0])
        nc_positions = np.array(positions)

        # define list of bonds
        start_index = 0  # index to start when assigning cage bonds
        bond_list = []
        if nc_at_origin:
            start_index = 1
            for i in range(1, num_cage_ncs+1):
                bond_list.append((0, i, 0))
                
        for i in range(start_index, start_index+(num_cage_ncs-1)):
            bond_list.append((i, i + 1, start_index))
        if num_cage_ncs > 2:
            bond_list.append((start_index+num_cage_ncs-1, start_index, start_index))

        # set ratios for analysis purposes
        bond_number_ratio = 1
        if num_cage_ncs == 2:
            bond_number_ratio = 1/2

        bond_length_ratio = 2*np.sin(np.pi/num_cage_ncs)

        super(PlanarConfigData, self).__init__(nc_positions, bond_list, bond_length_ratio, bond_number_ratio, num_cage_ncs)
