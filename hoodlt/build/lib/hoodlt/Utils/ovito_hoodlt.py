"""
:module: Visualize_files using the ovito software
:platform: Unix, Windows
:synopsis: Utility to visualize lattice structures

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, October2016
"""

import numpy as np
from hoodlt.Utils.LatticeCoords import LatCords


def make_ovito(lat):
    """Provides a file in extended xyz format for the lattice

    :param lat: :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
    :return: writes a new file in the current directory in extended xyz format that can be read by ovito
    """

    fn = lat.name()[0] + '.xyz'

    fid = open(fn, 'w')

    fid.write('%d \n' % np.sum(lat.typ))

    txt = 'Lattice="'

    ltc = LatCords(lat)

    for ind in range(3):
        vec = lat.get_a(ind)
        txt += '%1.7f %1.7f %1.7f ' % (vec[0], vec[1], vec[2])
    txt += '" Properties:pos:R:3:Radius:R:1'
    fid.write(txt + '\n')

    for ind_1 in range(len(lat.typ)):
        nm = chr(65+ind_1)
        rad = lat.radius[ind_1]
        for ind_2 in range(lat.typ[ind_1]):
            vec = ltc.fractional_unitcell(lat.get_v([ind_1, ind_2]))
            txt = nm + ' %1.7f %1.7f %1.7f %1.7f' % (vec[0], vec[1], vec[2], rad)
            fid.write(txt+'\n')

    fid.close()
