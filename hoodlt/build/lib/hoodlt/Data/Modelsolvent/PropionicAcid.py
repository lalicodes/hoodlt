"""
:module: PropionicAcid
:platform: Unix. Windows
:synopsis: Implements a class defining Propionic Acid Molecule

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> August 2019
.. history:
"""

from __future__ import division
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs
import numpy as np


class PropionicAcid(SolventAbs):
    """
    Defines solvent Propionic Acid (:math:`\\mbox{C}_3\\mbox{H}_6\\mbox{O}_2`)
    """

    def __init__(self, ff):
        """
        create the Propionic Acid molecule
        """

        super(PropionicAcid, self).__init__(1, ff, 11, 'PropionicAcid')  # always only 1 repeat

        self.mass[0] = self.ff_reader.get_molecular_weight('C136')
        self.mass[1] = self.ff_reader.get_molecular_weight('C267')
        self.mass[2] = self.ff_reader.get_molecular_weight('O269')
        self.mass[3] = self.ff_reader.get_molecular_weight('O268')
        self.mass[4] = self.ff_reader.get_molecular_weight('H270')
        self.mass[5] = self.ff_reader.get_molecular_weight('C135')
        self.mass[6:11] = self.ff_reader.get_molecular_weight('H140')

        self.charge[0] = self.ff_reader.get_charge('C136')
        self.charge[1] = self.ff_reader.get_charge('C267')
        self.charge[2] = self.ff_reader.get_charge('O269')
        self.charge[3] = self.ff_reader.get_charge('O268')
        self.charge[4] = self.ff_reader.get_charge('H270')
        self.charge[5] = self.ff_reader.get_charge('C135')
        self.charge[6:11] = self.ff_reader.get_charge('H140')

        # types
        self.types = ['C136', 'C267', 'O269', 'O268', 'H270', 'C135', 'H140']
        self.typeid = ['C136', 'C267', 'O269', 'O268', 'H270', 'C135'] + ['H140'] * 5

        # particle positions, in Angstroms
        self.position = self.build_chain()

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[5]  # points from one Carbon to another

    def build_chain(self):
        """
        called by initializer to get positions

        :return: position array for the Propionic Acid molecule
        """
        chain = np.zeros((11, 3))

        r_cc1 = self.ff_reader.get_bond_r0('C267-C136')
        r_oc1 = self.ff_reader.get_bond_r0('O269-C267')
        r_oc2 = self.ff_reader.get_bond_r0('O268-C267')
        r_ho = self.ff_reader.get_bond_r0('H270-O268')
        r_cc2 = self.ff_reader.get_bond_r0('C135-C136')
        r_hc1 = self.ff_reader.get_bond_r0('H140-C136')
        r_hc2 = self.ff_reader.get_bond_r0('H140-C135')

        theta_ccc = self.ff_reader.get_angle_t0('C267-C136-C135')
        theta_occ1 = self.ff_reader.get_angle_t0('O269-C267-C135')
        theta_coh = self.ff_reader.get_angle_t0('C267-O268-H270')
        theta_cch = self.ff_reader.get_angle_t0('C267-C136-H140')
        theta_occ2 = self.ff_reader.get_angle_t0('C136-C267-O268')
        theta_hch1 = self.ff_reader.get_angle_t0('H140-C136-H140')
        theta_hch2 = self.ff_reader.get_angle_t0('H140-C135-H140')
        theta_cch2 = self.ff_reader.get_angle_t0('C136-C135-H140')

        # equilibrium dihedral angle

        # supplement angles
        phi0 = np.pi - theta_ccc
        phi1 = np.pi - theta_occ1
        phi2 = np.pi - theta_occ2

        chain[1] = [r_cc1, 0, 0]
        chain[2] = [-r_oc1*np.cos(theta_occ1) + r_cc1, r_oc1*np.sin(phi1), 0]
        chain[3] = [-r_oc2*np.cos(theta_occ2) + r_cc1, -r_oc2*np.sin(phi2), 0]
        chain[4] = [r_cc1 + r_oc2*np.cos(np.pi-theta_occ2) + r_ho*np.cos(theta_coh-np.pi/2 - (np.pi/2-(np.pi-theta_occ2))),
                    -r_oc2*np.sin(np.pi-theta_occ2) - r_ho*np.sin(theta_coh-np.pi/2 - (np.pi/2-(np.pi-theta_occ2))), 0]
        chain[5] = [-r_cc2*np.cos(phi0), r_cc2*np.sin(phi0), 0]

        chain[6] = [-r_hc1*np.cos(theta_hch1/2)*np.sin(np.arccos(np.cos(theta_cch)/np.cos(theta_hch1/2))-np.pi/2),
                    -r_hc1*np.cos(theta_hch1/2)*np.cos(np.arccos(np.cos(theta_cch)/np.cos(theta_hch1/2))-np.pi/2),
                    r_hc1*np.sin(theta_hch1/2)]
        chain[7] = [-r_hc1*np.cos(theta_hch1/2)*np.sin(np.arccos(np.cos(theta_cch)/np.cos(theta_hch1/2))-np.pi/2),
                    -r_hc1*np.cos(theta_hch1/2)*np.cos(np.arccos(np.cos(theta_cch)/np.cos(theta_hch1/2))-np.pi/2),
                    -r_hc1*np.sin(theta_hch1/2)]

        ang0 = theta_cch2 - np.pi/2. - (np.pi/2.-phi0)
        ang1 = np.arccos(np.cos(theta_hch2)/np.cos(theta_hch2/2.))

        chain[8] = [-r_cc2*np.cos(phi0) - r_hc2*np.cos(ang0), r_cc2*np.sin(phi0) + r_hc2*np.sin(ang0), 0]
        chain[9] = [r_hc2*np.cos(theta_hch2/2.)*np.cos(np.pi-ang1-ang0) - r_cc2*np.cos(phi0),
                    r_cc2*np.sin(phi0) + r_hc2*np.cos(theta_hch2/2.)*np.sin(np.pi-ang1-ang0),
                    r_hc2*np.sin(theta_hch2/2.)]
        chain[10] = [r_hc2*np.cos(theta_hch2/2.)*np.cos(np.pi-ang1-ang0) - r_cc2*np.cos(phi0),
                     r_cc2*np.sin(phi0) + r_hc2*np.cos(theta_hch2/2.)*np.sin(np.pi-ang1-ang0),
                     -r_hc2*np.sin(theta_hch2/2.)]

        return chain

    def get_bonds(self):
        """
        override of method in BasicSystemEntity
        returns a list of the particle indices (2 particles per bond) involved in the bonds in the system

        :return: a 5x2 numpy array
        """

        return np.array([[1, 0], [2, 1], [3, 1], [4, 3], [5, 0], [6, 0], [7, 0], [8, 5], [9, 5], [10, 5]])

    def get_angles(self):
        """
        override of method in BasicSystemEntity
        returns a list of particle indices (3 particles per angle) involved in the angles in the system

        :return: a 6x3 numpy array
        """

        return np.array([[0, 1, 2], [0, 1, 3], [1, 3, 4], [1, 0, 5], [1, 0, 6], [1, 0, 7], [0, 5, 8], [0, 5, 9],
                         [0, 5, 10], [9, 5, 10], [6, 0, 7], [5, 0, 6], [8, 5, 9], [2, 1, 3], [5, 0, 7], [8, 5, 10]])

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

        return np.array([[3, 1, 0, 2]])

    def get_bond_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the bond types in the system
        """

        return ['C267-C136', 'O269-C267', 'O268-C267', 'H270-O268', 'C135-C136'] + ['H140-C136'] * 2 + ['H140-C135'] * 3

    def get_angle_types(self):
        """
        override of the method in BasicSystemEntity

        :return: list of all the angle types in the system
        """

        return ['C136-C267-O269', 'C136-C267-O268', 'C267-O268-H270', 'C267-C136-C135'] + ['C267-C136-H140'] * 2 + \
               ['C136-C135-H140'] * 3 + ['H140-C135-H140', 'H140-C136-H140', 'C135-C136-H140', 'H140-C135-H140'] + \
               ['O269-C267-O268', 'C135-C136-H140', 'H140-C135-H140']

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

        return ['O268-C267-C136-O269']
