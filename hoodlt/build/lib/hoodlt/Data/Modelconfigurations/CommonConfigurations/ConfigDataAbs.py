"""
:module: ConfigDataAbs
:platform: Unix, Windows
:synopsis: Stores info to easily build complex configurations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np


class ConfigDataAbs:
    """
    Abstract class which stores all the relevant information for building complex configurations of ncs
    """

    def __init__(self, nc_positions, bond_list):
        """

        :param nc_positions: list of positions of the center of masses of the ncs
        :param bond_list: list of triples (i, j, k) which represent nc bond indexes of type index k
        """

        self.nc_positions = nc_positions
        self.bond_list = bond_list
