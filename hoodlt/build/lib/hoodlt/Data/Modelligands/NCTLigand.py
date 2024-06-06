"""
:module: NCTLigand
:platform: Unix, Windows
:synopsis: Implements the polymer version of the abstract ligand class for use in NanoComposite Tecton Simulations

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> September 2018
.. history::
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class NCTLigand(LigandAbs):
    """Defines a simple polymer Ligand with a Nanocomposite tecton (NCT) endgroup
    """

    def __init__(self, repeats, ff, nct_type='A'):
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
        self.mass[ind+1] = self.ff_reader.get_molecular_weight(nct_type)
        self.mass[ind+2:ind+5] = self.ff_reader.get_molecular_weight('NCTfl')

        # set rigid body stuff
        for i in range(1, self.num_particles):
            if i >= self.num_particles - 5:
                self.body[i] = self.num_particles - 5
                if i == self.num_particles - 5:
                    self.moment_inertia[i] = np.array([0, .54, .54]) # np.array([0, .32, .32])

        # types
        self.types = ['Gr', 'PS', nct_type, 'NCTfl']
        self.typeid = ['Gr'] + ['PS']*repeats + [nct_type] + ['NCTfl']*4
        self.repeats = repeats

        self.position = self.build_chain()

    def build_chain(self):
        """builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """
        chain = np.zeros((self.num_particles, 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polymer')
        nct_bond_length = self.ff_reader.get_bond_r0('nct_tether')

        current_pt = np.array([0.0, 0.0, 0.0])

        chain[0] = np.array([-1.0, 0.0, 0.0])*poly_bond_length + np.array([1.0, 0.0, 0.0])*poly_bond_length

        for i in range(1, self.repeats + 1):
            ind = i
            pt = np.array([1.0, 0.0, 0.0])*float(poly_bond_length)
            current_pt += pt
            chain[ind] = current_pt
        ind += 1
        pt = np.array([1.0, 0.0, 0.0])*float(nct_bond_length)
        current_pt += pt
        chain[ind] = current_pt

        chain[ind + 1] = current_pt + np.array([0.0, 1.0, 0.0])*0.5  #0.4
        chain[ind + 2] = current_pt + np.array([0.0, -1.0, 0.0])*0.5  #0.4
        chain[ind + 3] = current_pt + np.array([0.0, 0.0, 1.0])*0.5  #0.4
        chain[ind + 4] = current_pt + np.array([0.0, 0.0, -1.0])*0.5  #0.4

        sub_pt = np.array([0, 0, 0])
        for i in range(self.num_particles):
            chain[i] = np.subtract(chain[i], sub_pt)

        return chain

    def get_bond_types(self):
        """utilized to dump_gsd

        :return: list of all the bond types in the system
        """
        return ['polymer'] * (self.repeats) + ['nct_tether']

    def get_bonds(self):
        """overrides method in LigandAbs to dump GSD files

        :return: n by 2 list of indices that share bond potentials
        """

        bonds = np.zeros((self.repeats+1, 2))

        for i in range(0, self.repeats+1):
            bonds[i] = [i, i+1]

        return bonds

    def get_angles(self):
        return []

    def get_dihedrals(self):
        return []

