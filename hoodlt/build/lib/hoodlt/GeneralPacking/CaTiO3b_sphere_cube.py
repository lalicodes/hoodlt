"""
:module: CaTiO3b_sphere_cube
:platform: Unix, Windows
:synopsis: Defines the classes implementing a correction to the packing fraction for A:spheres and B: cubes CaTiO3

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np
import types as ty

import hoodlt.Lattices.Lat2.CaTiO3b_lat as lt
from hoodlt.GeneralPacking.pf_sphere_cube import pf_s_c


def catio3_sphere_cube(l_value, a_nn_e, gam):
    """Computes the packing fraction for sphere cube CaTiO3b lattice

    :param l_value: lattice size
    :param a_nn_e: spacing
    :param gam: :math:`\\gamma` value
    :returns: A lattice object with the additional member function pf_s_c calculating the pf for the sphere-cube
    """

    a_lat = 1.0

    # the lattice is the same as the standard NaCl
    lat = lt.LatCaTiO3Base5(l_value, a_nn_e, gam)

    if gam > lat.gamma_crit[1]:
        a_lat = gam / (np.sqrt(2) - 1)
    lat.redefine_lattice_constant(a_lat)

    lat.pf_s_c = ty.MethodType(pf_s_c, lat)

    lat.gamma_crit_s_c = np.array([0.0, lat.gamma_crit[1], 1])

    return lat
