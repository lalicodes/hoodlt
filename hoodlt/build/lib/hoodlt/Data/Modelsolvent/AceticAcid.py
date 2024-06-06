"""
:module: AceticAcid
:platform: Unix. Windows
:synopsis: Implements a class defining Acetic Acid Molecule

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> August 2019
.. history:
"""

from __future__ import division
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs
import numpy as np


class AceticAcid(SolventAbs):
    """
    Defines solvent Acetic Acid (:math:`\\mbox{C}_2\\mbox{H}_4\\mbox{O}_2`)
    """

    def __init__(self, ff):
        """
        create the Acetic Acid molecule
        """

        super(AceticAcid, self).__init__(1, ff, 8, 'AceticAcid')  # always only 1 repeat

        self.mass[0] = self.ff_reader.get_molecular_weight('C267')
        self.mass[1] = self.ff_reader.get_molecular_weight('O269')
        self.mass[2] = self.ff_reader.get_molecular_weight('O268')
        self.mass[3] = self.ff_reader.get_molecular_weight('C135')
        self.mass[4] = self.ff_reader.get_molecular_weight('H270')
        self.mass[5:8] = self.ff_reader.get_molecular_weight('H140')

        self.charge[0] = self.ff_reader.get_charge('C267')
        self.charge[1] = self.ff_reader.get_charge('O269')
        self.charge[2] = self.ff_reader.get_charge('O268')
        self.charge[3] = self.ff_reader.get_charge('C135')
        self.charge[4] = self.ff_reader.get_charge('H270')
        self.charge[5:8] = self.ff_reader.get_charge('H140')

        # types
        self.types = ['C267', 'O269', 'O268', 'C135', 'H270', 'H140']
        self.typeid = ['C267', 'O269', 'O268', 'C135', 'H270'] + ['H140'] * 3

        # particle positions, in Angstroms
        self.position = self.build_chain()

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[0] - self.position[3]  # points from one Carbon to another

    def build_chain(self):
        """
        called by initializer to get positions
        :return: position array for the Acetic Acid molecule
        """
        chain = np.zeros((8, 3))

        r_oc1 = self.ff_reader.get_bond_r0('O269-C267')
        r_oc2 = self.ff_reader.get_bond_r0('O268-C267')
        r_cc = self.ff_reader.get_bond_r0('C135-C267')
        r_ho = self.ff_reader.get_bond_r0('H270-O268')
        r_hc = self.ff_reader.get_bond_r0('H140-C135')

        theta_oco = self.ff_reader.get_angle_t0('O269-C267-O268')
        theta_occ1 = self.ff_reader.get_angle_t0('O269-C267-C135')
        theta_coh = self.ff_reader.get_angle_t0('C267-O268-H270')
        theta_cch = self.ff_reader.get_angle_t0('C267-C135-H140')
        theta_occ2 = self.ff_reader.get_angle_t0('O268-C267-C135')
        theta_hch = self.ff_reader.get_angle_t0('H140-C135-H140')

        # equilibrium dihedral angle
        theta_d = np.pi

        # supplement angles
        phi1 = np.pi - theta_occ1
        phi2 = np.pi - theta_occ2
        phi3 = np.pi - phi2 - np.pi/2

        chain[1] = [-r_oc1*np.cos(theta_occ1), r_oc1*np.sin(phi1), 0]
        chain[2] = [-r_oc2*np.cos(theta_occ2), -r_oc2*np.sin(phi2), 0]
        chain[3] = [-r_cc, 0, 0]
        chain[4] = [-r_oc2*np.cos(theta_coh)+r_ho, -r_ho*np.cos(theta_coh-phi3) - r_oc2*np.sin(phi2), 0]
        chain[5] = [r_hc*np.cos(theta_cch)-r_cc, -r_hc*np.cos(theta_cch-np.pi/2), 0]
        chain[6] = [r_hc*np.cos(theta_cch)-r_cc, r_hc*(1/np.sqrt(3))*np.sin(theta_hch/2), r_hc*(2/np.sqrt(3))*np.sin(theta_hch/2)]
        chain[7] = [r_hc*np.cos(theta_cch)-r_cc, r_hc*(1/np.sqrt(3))*np.sin(theta_hch/2), -r_hc*(2/np.sqrt(3))*np.sin(theta_hch/2)]

        return chain

    def get_bonds(self):
        """
        override of method in BasicSystemEntity
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system

        :return: a 5x2 numpy array
        """

        return np.array([[0, 1], [0, 2], [0, 3], [2, 4], [3, 5], [3, 6], [3, 7]])

    def get_angles(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (3 particles per angle) involved in the angles in the system

        :return: a 6x3 numpy array
        """

        return np.array([[1, 0, 2], [1, 0, 3], [0, 2, 4], [0, 3, 5], [0, 3, 6], [0, 3, 7], [2, 0, 3], [5, 3, 6], [5, 3, 7], [6, 3, 7]])

    def get_dihedrals(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (4 particles per dihedral) involved in the dihedrals in the system

        :return: a nx4 numpy array
        """

        return np.array([[4, 2, 0, 1], [4, 2, 0, 3]])

    def get_impropers(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (4 particles per improper) involved in impropers

        :return: nx4 numpy array
        """

        return np.array([[0, 1, 2, 3]])

    def get_bond_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the bond types in the system
        """

        return ['O269-C267', 'O268-C267', 'C135-C267', 'H270-O268']+ ['H140-C135'] * 3

    def get_angle_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the angle types in the system
        """

        return ['O269-C267-O268', 'O269-C267-C135', 'C267-O268-H270'] + ['C267-C135-H140'] * 3 + ['O268-C267-C135'] + ['H140-C135-H140'] * 3

    def get_dihedral_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the dihedral types in the system
        """

        return ['H270-O268-C267-O269', 'H270-O268-C267-C135']

    def get_improper_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all improper types in the system
        """

        return ['C267-O269-O268-C135']
