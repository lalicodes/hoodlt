"""
:module: NCTLigand
:platform: Unix, Windows
:synopsis: Implements the polymer version of the abstract ligand class for use in NanoComposite Tecton Simulations

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> September 2018
.. history::
..              Jianshe Xia <xiajs6075@iccas.ac.cn> Oct 2020
..              - Added the method for computing the moment of inertia
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class NCTLigand(LigandAbs):
    """Defines a simple polymer Ligand with a Nanocomposite tecton (NCT) endgroup
    """

    def __init__(self, repeats, ff, nct_type='cA'):
        """The constructor

        :param repeats: number of repeatable polymer units in chain
        :param ff: forcefield
        :param nct_typ: NCT type (String)
        """

        super(NCTLigand, self).__init__(repeats, ff, repeats+2+4, 'NCTLigand')

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('Gr')
        for i in range(repeats):
            ind = i + 1
            self.mass[ind] = self.ff_reader.get_molecular_weight('PS')
        self.mass[repeats+1] = self.ff_reader.get_molecular_weight(nct_type)
        self.mass[repeats+2:repeats+6] = self.ff_reader.get_molecular_weight('NCTfl')

        # types
        self.types = ['Gr', 'PS', nct_type, 'NCTfl']
        self.typeid = ['Gr'] + ['PS']*repeats + [nct_type] + ['NCTfl']*4
        self.repeats = repeats

        self.position = self.build_chain()

        self.rigid_mass = self.mass[repeats + 1:].reshape((5,1))
        self.rigid_position = self.position[repeats + 1:]
        self.rigid_shift(-1 * self.rigid_center_of_mass())
        self.temp_moment_inertia = self.rigid_moment_of_inertia() 

        # set rigid body stuff
        for i in range(1, self.num_particles):
            if i >= self.num_particles - 5:
                self.body[i] = self.num_particles - 5
                if i == self.num_particles - 5:
                    # np.array([0, .54, .54]) # np.array([0, .32, .32])
                    self.moment_inertia[i] = self.temp_moment_inertia

    def build_chain(self):
        """builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """
        chain = np.zeros((self.num_particles, 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polymer')
        nct_bond_length = self.ff_reader.get_bond_r0('nct_tether')
        poly_angle = self.ff_reader.get_angle_t0('PS-PS-PS')
        nct_angle = self.ff_reader.get_angle_t0('nct_angle')
        sigma = self.ff_reader.get_nonbonded_sigma('PS', 'PS')

        stvec1 = np.array([poly_bond_length * np.sin(poly_angle / 2.0), poly_bond_length * np.cos(poly_angle / 2.0), 0.0])
        stvec2 = np.array([poly_bond_length * np.sin(poly_angle / 2.0), -poly_bond_length * np.cos(poly_angle / 2.0), 0.0])
        gr_pos = np.array([-poly_bond_length * np.sin(poly_angle / 2.0), poly_bond_length * np.cos(poly_angle / 2.0), 0.0])

        current = np.array([0.0, 0.0, 0.0])
        chain[0] = gr_pos
        chain[1] = current

        for i in range(2, self.repeats + 1):
            if i % 2 == 0:
                current = np.add(current, stvec1)
            else:
                current = np.add(current, stvec2)

            chain[i] = current

        pt = np.array([nct_bond_length * np.sin(nct_angle / 2.0), nct_bond_length * np.cos(nct_angle / 2.0), 0.0])
        current_pt = current + pt
        chain[self.repeats + 1] = current_pt
        ind = self.repeats + 1
        chain[ind + 1] = current_pt + np.array([0.0, 1.0, 0.0]) * sigma/2.0
        chain[ind + 2] = current_pt + np.array([0.0, -1.0, 0.0]) * sigma/2.0
        chain[ind + 3] = current_pt + np.array([0.0, 0.0, 1.0]) * sigma/2.0
        chain[ind + 4] = current_pt + np.array([0.0, 0.0, -1.0]) * sigma/2.0

        for i in range(self.num_particles):
            chain[i] = np.subtract(chain[i], gr_pos)

        return chain

    def get_bond_types(self):
        """utilized to dump_gsd

        :return: list of all the bond types in the system
        """
        return ['polymer'] * (self.repeats) + ['nct_tether']

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain 
        """
        return ['PS-PS-PS'] * (self.repeats - 1) + ['nct_angle']

    def get_bonds(self):
        """overrides method in LigandAbs to dump GSD files

        :return: n by 2 list of indices that share bond potentials
        """

        bonds = np.zeros((self.repeats+1, 2))

        for i in range(0, self.repeats+1):
            bonds[i] = [i, i+1]

        return bonds

    def get_angles(self):
        """
        used by other classes to add angles when the gsd is dumped
        :return: an n by 3 list of indexes that have an angle together	
        """

        angles = np.zeros((self.repeats, 3))
        for i in range(0, self.repeats):
            angles[i] = [i, i+1, i+2]

        return angles

    def get_dihedrals(self):
        return []

    def rigid_moment_of_inertia(self):
        """returns the moment of inertia
        :param: ctol tolerance for center of mass
        :param: dtol tolerance for nonzero diagonals
        :return: moment of inertia matrix (3,3)
        """
        i_diag = np.sum(self.rigid_mass * self.rigid_position ** 2)
        tensor = i_diag * np.eye(3, 3) - np.tensordot(self.rigid_mass * self.rigid_position, self.rigid_position,
                                                      axes=[0, 0])
        return np.array([tensor[0, 0], tensor[1, 1], tensor[2, 2]])

    def rigid_shift(self, vector):
        """adds the vector to all the positions in the solvent molecule
        :return: void
        """
        self.rigid_position = np.add(self.rigid_position, vector)

    def rigid_center_of_mass(self):
        """
        Calculate the center of mass of the rigid end. This function will be used in the constructor to ensure that the center
        of mass is located at the origin.
        :return: a numpy array containing the position of the center of mass
        """
        return np.sum(self.rigid_position * self.rigid_mass, axis=0) / np.sum(self.rigid_mass)
