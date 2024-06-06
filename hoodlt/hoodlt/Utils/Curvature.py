"""
:module: Curvature
:platform: Unix, Windows
:synopsis: Utility to calculate the Gaussian curvature of a polyhedra

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2017
"""

import numpy as np
import numpy.linalg as la


class CurvPolyHedra(object):
    """
    Class used to compute the curvature of a polyhedra
    
    :ivar qdih: dihedral angle
    :ivar coxeter_hc: Coxeter value
    :ivar c_topo: curvature calculated using topological formula (assumes all edges have the same length)
    :ivar c_metric: curvature calculated using the metric formula. 
    :ivar ico: number of icosahedral sites
    :ivar octa: number of 2pi/3 disclination sides
    :ivar edge length: average edge length
    """

    def __init__(self, poly):
        """The constructor

        :param poly: :class:`hoodlt.Utils.GeomClasses.PolyHedra`
        """

        # dihedral angle
        self.qdih = np.arccos(1.0/3.0)

        # coeff
        coeff = np.array([2*np.pi, -self.qdih])

        # Coxeter honeycomb value
        self.coxeter_hc = - coeff[0]/coeff[1]

        # curvature calculated using topological formula (assumes all edges have the same length)
        self.c_topo = np.zeros([2, 2])
        # curvature calculated using metric formula (edges have different lengths)
        self.c_metric = np.zeros([2, 2])
        # first index denotes the term with coefficients of pi and the one with coefficients of arccos(1/3)
        # the second index indicates 2pi/5q either 2pi/3 disclinations

        # reference curvature (curvature if all the lengths were the same)
        self.curv_reference_topo = (2 * np.pi - 5 * self.qdih) * poly.num_faces
        self.curv_reference_metric = (2 * np.pi - 5 * self.qdih) * poly.num_faces

        # number of icosahedral sites
        self.ico = 0

        # number of 2pi/3 disclination sides
        self.octa = 0

        for fc in poly.faces:
            n_edges = len(fc)
            cmass = np.sum(poly.vertices[fc], axis=0) / float(n_edges)

            self.c_topo[0, 0] += 1
            self.c_topo[1, 0] += n_edges
            self.c_metric[0, 0] += la.norm(cmass)
            self.c_metric[1, 0] += la.norm(cmass)*n_edges

            if n_edges == 5:
                self.ico += 1

        self.edge_length = self.c_metric[0, 0]/float(poly.num_faces)
        self.curv_reference_metric *= self.edge_length

        self.l_disp = 0.0
        for fc in poly.faces:
            n_edges = len(fc)
            cmass = np.sum(poly.vertices[fc], axis=0) / float(n_edges)
            ddist = la.norm(cmass)-self.edge_length
            self.l_disp += ddist**2/float(poly.num_faces)

        for ind_v in range(poly.coordination.shape[0]):
            if poly.coordination[ind_v] == 4:
                self.c_topo[0, 1] += -2
                self.c_topo[1, 1] += -8

                self.c_metric[0, 1] += -2 * self.edge_length
                self.c_metric[1, 1] += -8 * self.edge_length
                self.octa += 1

        self.curv_topo = np.sum(np.sum(self.c_topo, axis=1)*coeff)
        self.curv_metric = np.sum(np.sum(self.c_metric, axis=1)*coeff)
        # calculate dispersion of lengths (relative to average length)
        self.l_disp = np.sqrt(self.l_disp)/self.edge_length
