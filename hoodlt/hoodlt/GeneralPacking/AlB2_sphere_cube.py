"""
:module: AlB2_sphere_cube
:platform: Unix, Windows
:synopsis: Defines the classes implementing a correction to the packing fraction for A:spheres and B: cubes AlB2

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np
import types as ty

from hoodlt.GeneralPacking.pf_sphere_cube import pf_s_c

import hoodlt.Lattices.Lat2.AlB2_lat as lt


def alb2_sphere_cube(l_value, a_nn_e, gam):
    """Computes the packing fraction for sphere cube AlB2 lattice

    :param l_value: lattice size
    :param a_nn_e: spacing
    :param gam: :math:`\\gamma` value
    :returns: A lattice object with the additional member function pf_s_c calculating the pf for the sphere-cube
    """

    gam_1 = (4*np.sqrt(2)+np.sqrt(3) - np.sqrt(-1+8*np.sqrt(6)))/9.0
    gam_2 = 1/(2*np.sqrt(2))
    gam_3 = 2/(3*np.sqrt(3))
    gam_4 = 1/np.sqrt(3)

    lat = lt.LatAlB2Base3(l_value, a_nn_e, gam)

    if gam <= gam_1:
        a_lat = 1.0
        c_lat = 1.0
    elif gam_1 < gam < gam_2:
        # a_lat = (4*np.sqrt(2)*lat.gam+np.sqrt(3)*lat.gam+np.sqrt(21+(8*np.sqrt(6)-28)*lat.gam**2))/7.0
        a_lat = 1.0
        c_lat = gam/np.sqrt(3)+np.sqrt(-1+8*np.sqrt(2)*gam-8*gam**2)/np.sqrt(3)
    elif gam_2 <= gam < gam_3:
        a_lat = 2*np.sqrt(2)*gam
        c_lat = gam/np.sqrt(3)+np.sqrt(1-8*gam**2/3)
    elif gam_3 <= gam < gam_4:
        a_lat = 2*np.sqrt(2)*gam
        c_lat = 1.0
    else:
        a_lat = 2*np.sqrt(2)*gam
        c_lat = np.sqrt(3)*gam

    lat.redefine_lattice_constant([a_lat, a_lat, c_lat])

    lat.pf_s_c = ty.MethodType(pf_s_c, lat)

    lat.gamma_crit_s_c = np.array([0.0, gam_1, gam_2, gam_3, gam_4, 1])

    return lat


def alb2_sphere_cube_o2_sol1(l_value, a_nn_e, gam):
    """Computes the packing fraction for sphere cube AlB2 lattice orientation 2 solution 1

    :param l_value: lattice size
    :param a_nn_e: spacing
    :param gam: :math:`\\gamma` value
    :returns: A lattice object with the additional member function pf_s_c calculating the pf for the sphere-cube
    """

    gam_1 = (4*np.sqrt(2)+np.sqrt(3) - np.sqrt(-1+8*np.sqrt(6)))/9.0
    gam_2 = np.sqrt(2*(11+2*np.sqrt(6))/97)
    gam_3 = 1/np.sqrt(2)

    lat = lt.LatAlB2Base3(l_value, a_nn_e, gam)

    if gam <= gam_1:
        a_lat = 1.0
        c_lat = 1.0
    elif gam_1 < gam <= gam_2:
        a_lat = (4*np.sqrt(2)*lat.gam+np.sqrt(3)*lat.gam+np.sqrt(21+(8*np.sqrt(6)-28)*lat.gam**2))/7.0
        c_lat = a_lat
    elif gam_2 < gam <= gam_3:
        a_lat = 1.5*np.sqrt(2)*lat.gam
        c_lat = (lat.gam+np.sqrt(3-2*lat.gam**2))/np.sqrt(3)
    else:
        a_lat = 1.5*np.sqrt(2)*lat.gam
        c_lat = np.sqrt(3)*lat.gam

    lat.redefine_lattice_constant([a_lat, a_lat, c_lat])

    lat.pf_s_c = ty.MethodType(pf_s_c, lat)

    lat.gamma_crit_s_c = np.array([0.0, gam_1, gam_2, gam_3, 1])

    return lat


def alb2_sphere_cube_o2_sol2(l_value, a_nn_e, gam):
    """Computes the packing fraction for sphere cube AlB2 lattice orientation 2 solution 2

    :param l_value: lattice size
    :param a_nn_e: spacing
    :param gam: :math:`\\gamma` value
    :returns: A lattice object with the additional member function pf_s_c calculating the pf for the sphere-cube
    """

    gam_1 = (4*np.sqrt(2)+np.sqrt(3) - np.sqrt(-1+8*np.sqrt(6)))/9.0
    gam_2 = 1 / np.sqrt(3)
    gam_3 = 1/np.sqrt(2)

    lat = lt.LatAlB2Base3(l_value, a_nn_e, gam)

    if gam <= gam_1:
        a_lat = 1.0
        c_lat = 1.0
    elif gam_1 < gam <= gam_2:
        a_lat = 0.5*(2 * np.sqrt(2) * lat.gam + np.sqrt(2 * np.sqrt(3) * lat.gam - lat.gam ** 2))
        c_lat = 1.0
    elif gam_2 < gam <= gam_3:
        a_lat = np.sqrt(2)*lat.gam + np.sqrt(0.75-lat.gam**2)
        c_lat = np.sqrt(3)*lat.gam
    else:
        a_lat = 1.5*np.sqrt(2)*lat.gam
        c_lat = np.sqrt(3)*lat.gam

    lat.redefine_lattice_constant([a_lat, a_lat, c_lat])

    lat.pf_s_c = ty.MethodType(pf_s_c, lat)

    lat.gamma_crit_s_c = np.array([0.0, gam_1, gam_2, 1])

    return lat
