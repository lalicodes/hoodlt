"""
:module: NaCl_sphere_cube
:platform: Unix, Windows
:synopsis: Defines the classes implementing a correction to the packing fraction for A:spheres and B: cubes NaCl

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import types as ty
import numpy as np

from hoodlt.GeneralPacking.pf_sphere_cube import pf_s_c

import hoodlt.Lattices.Lat2.NaCl_lat as lt


def nacl_sphere_cube(l_value, a_nn_e, gam):
    """Computes the packing fraction for sphere cube NaCl lattices

    :param l_value: lattice size
    :param a_nn_e: spacing
    :param gam: :math:`\\gamma` value
    :returns: A lattice object with the additional member function pf_s_c calculating the pf for the sphere-cube
    """

    # the lattice is the same as the standard NaCl
    lat = lt.LatNaClBase8(l_value, a_nn_e, gam)

    lat.pf_s_c = ty.MethodType(pf_s_c, lat)

    lat.gamma_crit_s_c = np.array([0.0, lat.gamma_crit[1], 1])

    return lat
