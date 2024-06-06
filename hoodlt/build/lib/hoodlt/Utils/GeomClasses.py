"""
:module: GeomClasses
:platform: Unix, Windows
:synopsis: Utility to calculate different geometric properties of a lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2015
"""

import numpy as np


class PolyHedra(object):
    """
    Defines a simple class to parameterize a polyhedra
    """

    def __init__(self, vert, faces, qcor):
        """Constructor of the polyhedra

        :param vert: numpy array defining the vertices
        :param faces: list of integer numpy arrays defining the faces
        :param qcor: coordination of each vertex
        """

        self.vertices = vert
        self.faces = faces
        self.coordination = qcor
        self.num_faces = len(faces)
        self.num_vertices = vert.shape[0]
        self.num_edges = self.num_faces + self.num_vertices - 2


class VoroData(object):
    """
    Class used to store a Voronoi triangulation
    
    :ivar u_c: Unit cell object :class:`hoodlt.Utils.Lattice_UnitCell.UnitCell`
    :ivar num_vts: number of basis in the unit cell
    :ivar polyhedra_type: a numpy array containing # faces, # edges, # vertices for the given base
    :ivar vertices: a list of vectors defining the points of the voronoi cell for a given base
    :ivar regions: a list of arrays containing the indices defining each face for the cell for a given base
    :ivar neighs: a list of arrays containing the coordinates of the nearest neighbors of a given base
    :ivar neighs_rad: a list of arrays containing the radius of the nearest neighbors of a given base
    
    Examples
        .. code-block:: python
        
            from hoodlt.Utils.GeomClasses import VoroData
            from hoodlt.Utils.Lattice_UnitCell import UnitCell
            
            lat_pnt = UnitCell(lat)
            vdat = VoroData(lat_pnt)
            ind = 1
            ws = vdat.ws(ind)
            # returns the wigner seitz cell for the basis vector ind.
    """

    def __init__(self, u_c):
        """The constructor

        :param u_c: UnitCell object
        """

        self.u_c = u_c
        self.num_vts = self.u_c.num_basis
        self.points = np.zeros([self.num_vts, 3])

        ind_basis = 0
        for ind_1 in range(len(u_c.lat.typ)):
            for ind_2 in range(u_c.lat.typ[ind_1]):
                self.points[ind_basis] = self.u_c.pnts[self.u_c.basis_coord[ind_1][ind_2]]
                ind_basis +=1

        self.polyhedra_type = np.zeros([self.num_vts, 3], dtype=int)
        self.vertices = self.num_vts * [None]
        self.regions = self.num_vts * [None]
        self.neighs = self.num_vts * [None]
        self.neighs_id = self.num_vts *[None ]
        self.neighs_rad = self.num_vts * [None]

    def qcor(self, ind_vec, ind):
        """Returns the coordination of the vertex ind
        
        :param ind_vec: labels the corresponding vertex
        :param ind: basis index
        :return: coordination number as an integer
        :rtype: integer
        """

        qv = 0
        for reg in self.regions[ind]:
            qv += np.sum(np.equal(reg, ind_vec))

        return qv

    def ws(self, ind):
        """Returns the Wigner-Seitz cell as a PolyHedra object for basis ind

        :param ind: basis (integer)
        :return: The given polyhedra
        :rtype: PolyHedra
        """

        qcord = np.zeros(self.vertices[ind].shape[0])

        for ind_vec in range(self.vertices[ind].shape[0]):
            qcord[ind_vec] = self.qcor(ind_vec, ind)

        return PolyHedra(self.vertices[ind], self.regions[ind], qcord)

    def ws_simple(self, ind):
        """Returns the Wigner-Seitz cell as a list

            :param ind: basis (integer)
            :return: the given polyhedra
            :rtype: list
        """

        ind_total = self.u_c.u_cell[ind, 0]
        ind_1 = self.u_c.indices_all[ind_total, 0]
        rad = self.u_c.lat.radius[ind_1]

        return [self.points[ind], rad, self.vertices[ind], self.regions[ind]]
