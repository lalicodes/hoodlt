"""
:module: HydroCarbon
:platform: Unix, Windows
:synopsis: returns general hydrocarbon properties

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2016
"""

import numpy as np
import numpy.linalg as la


class HydroCarbon(object):
    """
    Defines the class defining the parameters of a hydrocarbon object
    """

    def __init__(self, num, unsat=None, amine=False):
        """
        Defines the constructor

        :param num: number of Carbons
        :param unsat: List containing saturation positions
        :param amine: True if the hydrocarbon terminates with an amine group

        ..note:: unsat is a list = [n1, n2, ..] n1 >= 1 is the position of the first double bond, etc..
        """

        if unsat is None:
            self.unsat = []
        else:
            self.unsat = unsat

        self.n = num
        # carbon end radius (in nm)
        self.c_end = 0.2
        # carbon single-bond distance (in nm)
        self.c_c = 0.153
        # carbon dihedral angle
        # self.c_dihed = 112.0
        self.c_dihed = 109.5
        self.c_angle = np.deg2rad(180-self.c_dihed)
        # length of a simple bond
        self.bond = self.c_c*np.cos(0.5*self.c_angle)
        # carbon double bond distance (in nm)
        self.bond_double = 0.1333
        self.d_angle = np.deg2rad(24.5)
        # maximum area density in inverse nm square
        self.area_dens_max = 5.5
        # length of amine terminated
        self.amine = self.c_end

        self.end = self.c_end

        if amine:
            self.end = self.amine

        self.n_double = len(self.unsat)
        # splice the chain into segments of saturated chains
        self.segs = (self.n_double + 1)*[0]

        if not self.unsat:
            self.segs[0] = [0, self.n]
        else:
            self.segs[0] = [0, self.unsat[0]]
            for ind in range(len(self.unsat)-1):
                self.segs[ind+1] = [self.unsat[ind] + 1, self.unsat[ind+1]]
            self.segs[-1] = [self.unsat[-1] + 1, self.n]

    def max_length(self):
        """
        Returns the maximum length of a hydrocarbon chain
        :return: hydrocarbon length (in nanometers)
        :rtype: double
        """

        sat_l = 0.0

        for ind, val in enumerate(self.segs):
            sat_l += val[1] - val[0]

        return self.bond*sat_l + len(self.unsat)*self.bond_double + self.end

    def unsat_length(self, n_c):
        """returns the length of a saturated length with num_c carbons (excluding the ends)

        :param n_c: number of carbons within the saturated chain
        :return: length
        :rtype: double
        """
        return n_c*self.bond

    def max_extent(self):
        """returns the maximum end to end distance

        :return: maximum end to end distance (in nanometers)
        :rtype: double
        """

        print("this result is only for saturated and oleic acid, other cases remain as TODO")

        if self.n_double == 0:
            e_to_e = self.max_length()
        elif self.n_double == 1:
            c_val = self.bond*np.sin(0.5*self.d_angle)*(self.segs[-1][1]-self.segs[-1][0])
            vec = np.array([((self.n-1)*self.bond+self.end)*np.cos(self.d_angle), c_val])
            e_to_e = la.norm(vec)
        else:
            raise ValueError('Saturations larger than 1 not computed yet')

        return e_to_e

    def lamb_matrix(self):
        """helper function to compute optimal distances as used in other applications

        :return: Returns a numpy array containing the maximum length and the maximum area density
        :rtype: numpy.ndarray
        """

        return np.array([self.max_length(), self.area_dens_max])
