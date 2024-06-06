"""
:module: AnalyzeBasicEntity
:Platform: Windows, Unix
:synopsis: Abstract class which contains methods to do calculations with individual basic system entities
nanoparticle

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, July 2019
.. history:
"""

import numpy as np
import numpy.linalg as la


class AnalyzeBasicEntity:
    """
    Abstract class which contains calculations that can be done on any basic system entity
    """

    def __init__(self, entity):
        """

        :param entity: BasicSystemEntity object
        """

        self.entity = entity

    def average_position(self):
        """
        Calculates the average position of the particles on the basic entity

        :return: the average position, as a numpy array of length 3
        """

        pos_avg = np.zeros(3)

        for pos in self.entity.position:
            pos_avg += pos

        return pos_avg / len(self.entity.position)

    def max_radius(self):
        """
        Calculates the maximum distance of any particle from the average position of the particles on the entity

        :return: radius of the entity
        """

        pos_avg = self.average_position()

        max_dist = 0
        for pos in self.entity.position:
            dist = la.norm(pos - pos_avg)
            if dist > max_dist:
                max_dist = dist

        return max_dist

    def position_by_type(self, type):
        """
        Creates position array by specific particle type

        :param type: string corresponding to specific particle type needed
        :return: numpy array containing all positions of requested type
        """

        pos_list = []
        for ind, typ in enumerate(self.entity.typeid):
            if typ == self.entity.types[self.entity.types.index(type)]:
                pos_list.append(self.entity.position[ind])

        return np.array(pos_list)
