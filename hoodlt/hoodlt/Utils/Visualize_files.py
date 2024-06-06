"""
:module: Visualize_files
:platform: Unix, Windows
:synopsis: Utility to visualize lattice structures

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, September2014
"""

import numpy as np
from numpy import linalg as la
import string


def make_hoomd(lat):
    """Provides a file of the lattice in hoomd form

    :param lat: DMixtLattice object
    :return: does not return anything but writes a new file in the current directory.
    """

    fn = lat.name()[0] + '.xml'

    print('preparing hoomd file')

    fid = open(fn, 'w')

    fid.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    fid.write('<hoomd_xml version="1.5">\n')
    fid.write('<configuration time_step="0" dimensions="3" natoms="%d">\n' % np.sum(lat.typ))

    lx = lat.a_axis
    ly = lat.b_axis
    lz = lat.c_axis
    [xy, xz, yz] = get_tilt_cell(lat.a_vector)

    print('tilt factors are %2.5f %2.5f %2.5f' % (xy, xz, yz))

    txt = '<box lx="%3.6f" ly="%3.6f" lz="%3.6f" xy="%3.6f" xz="%3.6f" yz="%3.6f"/>\n'
    fid.write(txt % (lx, ly, lz, xy, xz, yz))

    fid.write('<type>\n')
    for ind_l in range(len(lat.typ)):
        let = string.ascii_uppercase[ind_l]
        for ind_k in range(lat.typ[ind_l]):
            fid.write('%s\n' % let)
    fid.write('</type>\n')

    fid.write('<position>\n')
    for ind_l in range(len(lat.typ)):
        for ind_k in range(lat.typ[ind_l]):
            v = lat.get_v([ind_l, ind_k])
            fid.write('%2.5f %2.5f %2.5f\n' % (v[0], v[1], v[2]))
    fid.write('</position>\n')

    if lat.typ.size == 2:
        fid.write('<diameter>\n')
        for ind_l in range(len(lat.typ)):
            rad = 1.0
            if ind_l == 1:
                rad = lat.gam
            for ind_k in range(lat.typ[ind_l]):
                fid.write('%2.5f\n' % rad)
        fid.write('</diameter>\n')

    fid.write('</configuration>\n')
    fid.write('</hoomd_xml>\n')

    fid.close()


def make_pdb(lat):
    """Provides a file of the lattice in pdb form

    :param lat: DMixtLattice object
    :return: does not return anything but writes a new file in the current directory.
    """
    fn = lat.name()[0] + '.pdb'

    fid = open(fn, 'w')

    # write title
    fid.write('TITLE     %s %s\n' % (lat.name()[0], lat.name()[3]))

    # compute angles
    norm_a = la.norm(lat.get_a(0))
    norm_b = la.norm(lat.get_a(1))
    norm_c = la.norm(lat.get_a(2))

    alpha = np.rad2deg(np.arccos(np.dot(lat.get_a(0), lat.get_a(2))/(norm_a*norm_c)))
    beta = np.rad2deg(np.arccos(np.dot(lat.get_a(1), lat.get_a(2))/(norm_b*norm_c)))
    gamma = np.rad2deg(np.arccos(np.dot(lat.get_a(1), lat.get_a(0))/(norm_b*norm_a)))
    # write crystal structure
    txt = 'CRYST1    %2.1f   %2.1f   %2.1f  %2.1f  %2.1f  %2.1f %s\n'
    fid.write(txt % (norm_a, norm_b, norm_c, alpha, beta, gamma, lat.name()[2]))

    enm = 0
    for ind_l in range(len(lat.typ)):
        let = string.ascii_uppercase[ind_l]
        for ind_k in range(lat.typ[ind_l]):
            enm += 1
            vec = lat.get_v([ind_l, ind_k])
            fid.write('HETATM      %d  %s      %3.4f   %3.4f   %3.4f\n' % (enm, let, vec[0], vec[1], vec[2]))

    fid.close()


def get_tilt_cell(vec):
    """Returns the tilt vectors for a unit cell

    :param vec: a 3x3 numpy array containing the unit vectors
    :return: list of the three numbers as [xy, xz, yz]
    :rtype: list
    """
    lx = la.norm(vec[0])
    a_2x = np.dot(vec[0], vec[1])/lx
    ly = np.sqrt(la.norm(vec[1])**2-a_2x**2)
    xy = a_2x/ly

    vct = np.cross(vec[0], vec[1])
    lz = np.dot(vec[2], vct)/la.norm(vct)
    a_3x = np.dot(vec[0], vec[2])/lx
    xz = a_3x/lz
    yz = (np.dot(vec[1], vec[2])-a_2x*a_3x)/(lx*lz)

    return [xy, xz, yz]
