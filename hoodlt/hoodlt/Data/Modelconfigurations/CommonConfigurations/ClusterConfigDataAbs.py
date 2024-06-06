"""
:module: ClusterConfigDataAbs
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze a cluster configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

from hoodlt.Data.Modelconfigurations.CommonConfigurations.ConfigDataAbs import ConfigDataAbs


class ClusterConfigDataAbs(ConfigDataAbs):
    """
    Abstract class which stores all the relevant information for the cluster configurations
    """

    def __init__(self, nc_positions, bond_list, bond_length_ratio, bond_number_ratio, num_cage_ncs):
        """

        :param nc_positions: list of positions of the center of masses of the ncs
        :param bond_list: list of triples (i, j, k) which represent nc bonds between indexes i,j of type index k
        :param bond_length_ratio: ratio of the length of the cage-cage bonds to the length of the center-cage bonds
        :param bond_number_ratio: ratio of the number of cage-cage bonds to the number of center-cage bonds
        :param num_cage_ncs: the number of ncs in the cage
        """

        super(ClusterConfigDataAbs, self).__init__(nc_positions, bond_list)

        # constants for analysis, subclasses must override these fields
        self.bond_length_ratio = bond_length_ratio
        self.bond_number_ratio = bond_number_ratio
        self.num_cage_ncs = num_cage_ncs
