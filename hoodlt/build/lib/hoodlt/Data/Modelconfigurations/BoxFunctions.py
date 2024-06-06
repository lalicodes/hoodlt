"""
:module: BoxFunctions
:platform: Unix, Windows
:synopsis: class to manage box functions, implemented through freud.box
.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> May 2022
.. history:
"""

import numpy as np
import freud


class BoxFunctions(object):

    def __init__(self, l_b):
        """
        :param l_b: six dimensional box object
        creates several freud objectd
        """

        self.box = freud.box.Box(Lx=l_b[0], Ly=l_b[1], Lz=l_b[2], xy=l_b[3], xz=l_b[4], yz=l_b[5])

    def compute_distances(self, pos1, pos2):
        """
        return distances between pos1 and pos2

        :param pos1: numpy array of points
        :param pos2: numpy array of points
        :return : distances between points1 and points 2
        """

        return self.box.compute_distances(pos1, pos2)

    def compute_all_distances(self, pos1, pos2):
        """
        return distances between pos1 and pos2

        :param pos1: numpy array of points
        :param pos2: numpy array ofpoints
        :return : all distances between points1 and points 2
        """

        return self.box.compute_all_distances(pos1, pos2)

    def wrap(self, vec):
        """
        wraps the vectors within the box dimensions

        :param vec: vector array
        :return : wrapped vectors
        """

        return self.box.wrap(vec)

    def get_images(self, vec):
        """
        returns the images of vectors

        :param vec: vector
        """

        return self.box.get_images(vec)

    def unwrap(self, vec, imgs, out_val=None):
        """
        unwraps the vectors within the box dimensions

        :param vec: to-be-unwrapped vector array
        :param imgs: vector(s) image indices
        :param out_val: array with unwrapped vectors (if None, create an array)
        :return : unwrapped vectors
        """

        return self.box.unwrap(vec, imgs, out=out_val)

    def make_grid(self, tol):
        """
        return a list of point coordinates with tolerance

        :param: tol
        :return: numpy array with coordinates
        """

        vec_lat = np.array([self.box.Lx, self.box.Ly, self.box.Lz])
        vec_l = np.floor(vec_lat/np.array(tol)).astype(int)
        delta_l = 0.5/vec_l
        vec_m = [delta_l[ind]+np.arange(0, vec_l[ind])/vec_l[ind] for ind in range(3)]
        mgrid = [[ind1, ind2, ind3] for ind1 in vec_m[0] for ind2 in vec_m[1] for ind3 in vec_m[2]]
        return self.box.wrap(self.box.make_absolute(np.array(mgrid)))

    def lsize(self):
        """
        return box size
        """
        return self.box.L
