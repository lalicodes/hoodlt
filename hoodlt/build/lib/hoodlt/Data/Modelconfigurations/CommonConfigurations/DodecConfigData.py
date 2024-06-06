"""
:module: DodecConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build and analyze a dodecahedral configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ClusterConfigDataAbs import ClusterConfigDataAbs


class DodecConfigData(ClusterConfigDataAbs):
    """
    Class which stores all the data for building and analyzing a dodecahedral configuration of ncs
    """

    def __init__(self, dist_from_origin, nc_at_origin=True):
        """

        :param dist_from_origin: distance the cage ncs should be placed from the origin
        :param nc_at_origin: True if the position and bond list should include an nc at the origin
        """

        # define positions as scaled unit vectors
        positions = []

        # unit vectors
        d = 1 / np.sqrt(3)
        r = (1 + np.sqrt(5)) / 2
        n = 1 / 1.732  # normalization factor for vectors involving r
        unit_vectors = [(d,d,d),(d,d,-d),(d,-d,d),(d,-d,-d),(-d,d,d),(-d,d,-d),(-d,-d,d),(-d,-d,-d),
                        (0,n*r,n/r),(0,-n*r,n/r),(0,n*r,-n/r),(0,-n*r,-n/r),
                        (n/r,0,n*r),(n/r,0,-n*r),(-n/r,0,n*r),(-n/r,0,-n*r),
                        (n*r,n/r,0),(n*r,-n/r,0),(-n*r,n/r,0),(-n*r,-n/r,0)]

        if nc_at_origin:
            unit_vectors = [(0, 0, 0)] + unit_vectors

        for unit_vector in unit_vectors:
            positions.append([dist_from_origin*unit_vector[0],
                              dist_from_origin*unit_vector[1],
                              dist_from_origin*unit_vector[2]])

        nc_positions = np.array(positions)

        # define list of bonds
        num_cage_ncs = 20
        bond_list = []

        if nc_at_origin:
            for i in range(1, num_cage_ncs+1):
                bond_list.append((0, i, 0))

        # TODO add indexes for cage-cage bonds

        # set ratios for analysis purposes
        bond_number_ratio = 3/2
        bond_length_ratio = 4/np.sqrt(3)/(1+np.sqrt(5))

        super(DodecConfigData, self).__init__(nc_positions, bond_list, bond_length_ratio, bond_number_ratio, num_cage_ncs)
