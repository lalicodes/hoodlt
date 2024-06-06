"""
:module: Hydrazine
:platform: Unix. Windows
:synopsis: Implements a class defining Hydrazine

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> April 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Removed unecessary functions
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - edited get_vector() method
..                  - made class work with the new BasicSystemEntity class
"""

from __future__ import division
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs
import numpy as np

class Hydrazine(SolventAbs):
    """
    Defines solvent Hydrazine (:math:`\\mbox{N}_2\\mbox{H}_4`) 
    """

    def __init__(self, ff):
        """
        create the Hydrazine molecule
        """

        super(Hydrazine, self).__init__(1, ff, 6, 'Hydrazine')  # always only 1 repeat

        self.mass[0:2] = self.ff_reader.get_molecular_weight('NH2')
        self.mass[2:4] = self.ff_reader.get_molecular_weight('HA')
        self.mass[4:6] = self.ff_reader.get_molecular_weight('HB')

        self.charge[0:2] = self.ff_reader.get_charge('NH2')
        self.charge[2:4] = self.ff_reader.get_charge('HA')
        self.charge[4:6] = self.ff_reader.get_charge('HB')

        # types
        self.types = ['NH2', 'HA', 'HB']
        self.typeid = ['NH2']*2+['HA']*2+['HB']*2

        # particle positions, in Angstroms
        self.position = self.build_chain()

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]  # points from one NH2 to the other NH2

    def build_chain(self):
        """
        called by initializer to get positions

        :return: position array for the hydrocarbon chains
        """
        chain = np.zeros((6, 3))

        r_nn = self.ff_reader.get_bond_r0('NH2-NH2')
        r_na = self.ff_reader.get_bond_r0('NH2-HA')
        r_nb = self.ff_reader.get_bond_r0('NH2-HB')
        theta_ann = self.ff_reader.get_angle_t0('HA-NH2-NH2')
        theta_bnn = self.ff_reader.get_angle_t0('HB-NH2-NH2')
        theta_anb = self.ff_reader.get_angle_t0('HA-NH2-HB')

        #equilibrium dihedral angle
        theta_d = np.pi / 2
        #supplement of ann/bnn angles
        theta_bnn_sup = np.pi - theta_bnn
        theta_ann_sup = np.pi - theta_ann
        # some more angles
        phi_1 = theta_bnn - np.pi/2
        phi_2 = theta_ann - np.pi/2
        theta_p = np.arccos(np.cos(theta_anb)/np.cos(phi_2)*np.cos(phi_1))

        chain[1] = [0,0,r_nn]
        chain[2] = [r_na*np.sin(theta_ann), 0, r_na*np.cos(theta_ann)]
        chain[3] = [r_na*np.sin(theta_ann_sup)*np.cos(theta_d+theta_p), r_na*np.sin(theta_ann_sup)*np.sin(theta_d+theta_p), r_nn+r_na*np.cos(theta_ann_sup)]
        chain[4] = [r_nb*np.sin(theta_bnn)*np.cos(theta_p), r_nb*np.sin(theta_bnn)*np.sin(theta_p), r_nb*np.cos(theta_bnn)]
        chain[5] = [r_nb*np.sin(theta_bnn_sup)*np.cos(theta_d), r_nb*np.sin(theta_bnn_sup)*np.sin(theta_d), r_nn+r_nb*np.cos(theta_bnn_sup)]

        return chain

    def get_bonds(self):
        """
        override of method in BasicSystemEntity
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system
        
        :return: a 5x2 numpy array
        """

        return np.array([[0,1], [0,2], [1,3], [0,4], [1,5]])

    def get_angles(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (3 particles per angle) involved in the angles in the system
        
        :return: a 6x3 numpy array
        """

        return np.array([[2,0,1], [3,1,0], [1,0,4], [5,1,0], [3,1,5], [2,0,4]])

    def get_dihedrals(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (4 particles per dihedral) involved in the dihedrals in the system
        
        :return: a 4x4 numpy array
        """

        return np.array([[2,0,1,5], [3,1,0,4], [2,0,1,3], [4,0,1,5]])

    def get_bond_types(self):
        """
        override of the method in BasicSystemEntity
        
        :return: list of all the bond types in the system
        """

        return ['NH2-NH2'] + ['NH2-HA']*2 + ['NH2-HB']*2

    def get_angle_types(self):
        """
        override of the method in BasicSystemEntity
        
        :return: list of all the angle types in the system
        """

        return ['HA-NH2-NH2']*2 + ['HB-NH2-NH2']*2 + ['HA-NH2-HB']*2

    def get_dihedral_types(self):
        """
        override of the method in BasicSystemEntity
        
        :return: list of all the dihedral types in the system
        """

        return ['HA-NH2-NH2-HB']*2 + ['HA-NH2-NH2-HA'] + ['HB-NH2-NH2-HB']
