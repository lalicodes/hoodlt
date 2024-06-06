"""
:module: Yico_Group
:platform: Unix, Windows, OS
:synopsis: Defines The icosahedron group

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>  May 2017
"""

from __future__ import division
import numpy as np
import numpy.linalg as la
import hoodlt.Groups.SU as SU


class YGroup(object):
    """
    Defines the icosahedron group :math:`{\cal Y}` containing 120 elements
    """

    def __init__(self):
        """The Constructor
        """
        # group order
        self.order = 120

        # define edges
        edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 5), (1, 8), (1, 9), (1, 2), (2, 3), (2, 9),
                 (2, 10), (3, 4), (3, 10), (3, 6), (4, 5), (4, 6), (4, 7), (5, 7), (5, 8), (6, 7), (7, 8), (8, 9),
                 (9, 10), (10, 6), (6, 11), (7, 11), (8, 11), (9, 11), (10, 11)]

        # define faces
        faces = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1), (1, 5, 8), (1, 8, 9), (1, 9, 2), (2, 9, 10),
                 (2, 10, 3),
                 (3, 10, 6), (3, 6, 4), (4, 6, 7), (4, 7, 5), (5, 7, 8), (6, 10, 11), (6, 11, 7), (7, 11, 8),
                 (8, 11, 9), (9, 11, 10)]

        # define points (and axis)
        x = 1 / np.sqrt(5)
        y = 2.0 / np.sqrt(5)
        self.pnt_v = np.zeros([12, 3])
        self.pnt_v[0] = np.array([0.0, 0.0, 1.0])
        for ind in range(5):
            cn = np.cos(2 * ind * np.pi / 5.0)
            sn = np.sin(2 * ind * np.pi / 5.0)
            self.pnt_v[ind + 1] = np.array([y * cn, y * sn, x])
        self.pnt_v[6:11] = -self.pnt_v[1:6]
        self.pnt_v[-1] = -self.pnt_v[0]

        # define normal to edges
        self.pnt_e = np.zeros((30, 3))
        for ind, ed in enumerate(edges):
            self.pnt_e[ind] = self.pnt_v[ed[0]] + self.pnt_v[ed[1]]

        mod = np.zeros([30, 1])
        mod[:, 0] = la.norm(self.pnt_e, axis=1)
        self.pnt_e = self.pnt_e / mod

        # define normal to faces
        self.pnt_f = np.zeros((20, 3))
        for ind, pf in enumerate(faces):
            self.pnt_f[ind] = self.pnt_v[pf[0]] + self.pnt_v[pf[1]] + self.pnt_v[pf[2]]

        mod = np.zeros([20, 1])
        mod[:, 0] = la.norm(self.pnt_f, axis=1)
        self.pnt_f = self.pnt_f / mod

        # define the group classes and names
        self.classes = np.array([1, 1, 12, 12, 12, 12, 20, 20, 30], dtype=int)
        self.classes_cum = np.cumsum(self.classes)-1
        self.classes_name = ['C(0)', 'C(0,bar)', 'C(1,5)', 'C(2,5)', 'C(3,5)', 'C(4,5)', 'C(3,1)', 'C(3,2)', 'C(2)']

    def group_class(self, i, j):
        """Returns the element j of the class i
        
        :param i: class index
        :param j: element indes
        :return: Returns the element of the class
        :rtype: Quaternion
        """

        if i == 0:
            return SU.SU2(0.0, np.array([0.0, 0.0, 1.0])).rot

        elif i == 1:
            return SU.SU2(np.pi, np.array([0.0, 0.0, 1.0])).rot

        elif (i > 1) and (i < 6):
            ang = np.pi * (i - 1) / 5.0
            return SU.SU2(ang, self.pnt_v[j]).rot

        elif (i > 5) and (i < 8):
            ang = np.pi * (i - 5) / 3.0
            return SU.SU2(ang, self.pnt_f[j]).rot

        else:
            ang = 0.5*np.pi
            return SU.SU2(ang, self.pnt_e[j]).rot

    def flat_to_pair_index(self, i):
        """Returns the index of class from flat indexes from 0 to 119

                :param i: index element
                :return: a pair identifying class index
                :rtype: tuple
        """
        # map the element of the group to its corresponding class
        ind1 = int(np.digitize(i, self.classes_cum, right=True))
        if ind1 == 0:
            ind2 = 0
        else:
            ind2 = i - self.classes_cum[ind1 - 1] - 1

        return (ind1, ind2)

    def element_indx(self, i):
        """Returns group elements indexed from 0 to 119
        
        :param i: index element
        :return: The group element
        :rtype: Quaternion
        """

        ind1, ind2 = self.flat_to_pair_index(i)
        return self.group_class(ind1, ind2)
