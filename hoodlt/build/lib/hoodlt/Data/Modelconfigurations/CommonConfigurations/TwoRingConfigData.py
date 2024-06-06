"""
:module: TwoRingConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build a configuration of two rings of ncs

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ConfigDataAbs import ConfigDataAbs


class TwoRingConfigData(ConfigDataAbs):
    """
    Class which stores all the data necessary to build a configuration of ncs which has two rings of ncs centered around
    the z-axis, one above the x-y plane, and one below the x-y plane by the same amount
    """

    def __init__(self, num_ncs_per_ring, dist_from_origin):
        """

        :param num_ncs_per_ring: number of ncs to be spaced equally around in each ring
        :param dist_from_origin: distance the ncs should be from the origin
        """

        # set the positions of the ncs
        positions = []
        for i in range(2*num_ncs_per_ring):
            phi = np.pi / num_ncs_per_ring * i
            theta = np.pi / 2 - (np.power(-1, i)) * np.arctan(1 / 2)
            positions.append((dist_from_origin * np.sin(theta) * np.cos(phi),
                              dist_from_origin * np.sin(theta) * np.sin(phi),
                              dist_from_origin * np.cos(theta)))

        nc_positions = np.array(positions)

        # define a list of bonds
        bond_list = []
        for i in range(2*num_ncs_per_ring):
            bond_list.append((i, i+1, 0))

        bond_list.append((2*num_ncs_per_ring - 1, 0))

        super(TwoRingConfigData, self).__init__(nc_positions, bond_list)
