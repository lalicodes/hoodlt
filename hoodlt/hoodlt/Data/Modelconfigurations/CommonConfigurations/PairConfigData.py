"""
:module: PairConfigData
:platform: Unix, Windows
:synopsis: Stores info to easily build a pair configuration

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - bond list is now a triple with a bond type index as the 3rd entry
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added documentation
"""

import numpy as np
from hoodlt.Data.Modelconfigurations.CommonConfigurations.ConfigDataAbs import ConfigDataAbs

class PairConfigData(ConfigDataAbs):
    """
    Class that stores all the data for building a pair of ncs on the x-axis equidistant from the origin
    """

    def __init__(self, distance):
        """

        :param distance: The distance of each nc from the origin
        """

        bond_list = [(0, 1, 0)]
        nc_positions = np.array([[-distance/2, 0, 0], [distance/2, 0, 0]])

        super(PairConfigData, self).__init__(nc_positions, bond_list)
