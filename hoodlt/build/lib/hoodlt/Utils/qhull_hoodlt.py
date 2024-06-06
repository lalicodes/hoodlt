"""
:module: qhull_hoodlt
:platform: Unix, Windows
   :synopsis: Utility to calculate different geometric properties of a lattice using qhull python voronoi

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, December 2016
"""

import numpy as np

from scipy.spatial import Voronoi
from hoodlt.Utils.GeomClasses import VoroData
from hoodlt.Utils.Lattice_UnitCell import UnitCell

def qhull_hoodlt(lat):
    """Computes the Voronoi cell of a lattice

    :param lat: :class:`hoodlt.D_matrix.Dmatrix_Mixt_Lattices.DMixtLattice`
    :return: Voronoi parameters of the unit cell
    :rtype: VoroData
    """

    # this is the file that qhull reads
    lat_pnt = UnitCell(lat)

    vor = Voronoi(lat_pnt.pnts)
    vdat = VoroData(lat_pnt)

    ind_basis = 0
    for ind_1 in range(len(lat_pnt.lat.typ)):
        for ind_2 in range(lat_pnt.lat.typ[ind_1]):
            ip = lat_pnt.basis_coord[ind_1][ind_2]
            pnt_reg = vor.regions[vor.point_region[ip]]
            rdg = vor.ridge_vertices

            num_v = len(pnt_reg)
            ws = []
            vdat.vertices[ind_basis] = vor.vertices[pnt_reg]
            for v in rdg:
                if set(v).issubset(set(pnt_reg)):
                    # convert these points to relative coordinates
                    w = np.zeros([len(v)], dtype=int)
                    for ind, val in enumerate(v):
                        w[ind] = pnt_reg.index(val)
                    ws.append(w)

            num_f = len(ws)
            num_e = num_f + num_v - 2
            vdat.polyhedra_type[ind_basis] = np.array([num_f, num_e, num_v], dtype=int)

            vdat.regions[ind_basis] = ws
            ind_basis += 1

    return vdat
