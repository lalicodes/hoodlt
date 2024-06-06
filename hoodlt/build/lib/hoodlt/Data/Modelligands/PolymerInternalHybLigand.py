"""
:module: PolymerInternalHybLigand
:platform: Unix, Windows
:synopsis: Implements the polymer version of the abstract ligand class for use in  hydrogen-bonding simulations

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> September 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class PolymerInternalHybLigand(LigandAbs):
    """Defines the Polymer Ligand with hybridizable flanking sites
    """

    def __init__(self, repeats, ff):
        """The constructor

        :param repeats: number of repeatable polymer units in chain
        """

        super(PolymerInternalHybLigand, self).__init__(repeats, ff, 2*repeats+1, 'PolymerInternalHyb')

        # masses
        self.mass[0] = self.ff_reader.get_molecular_weight('Gr')
        for i in range(repeats):
            ind = i*2 + 1
            self.mass[ind] = self.ff_reader.get_molecular_weight('cPEG')
            self.mass[ind + 1] = self.ff_reader.get_molecular_weight('PEG')

        # types
        self.types = ['Gr', 'cPEG', 'PEG']
        self.typeid = ['Gr'] + ['cPEG', 'PEG']*repeats

        # set rigid body stuff
        for i in range(1, self.num_particles):
            if i%2 == 0:
                self.body[i] = i-1
            else:
                self.body[i] = i
                self.moment_inertia[i] = np.array([0, 0.625, 0.625])

        self.position = self.build_chain()

    def build_chain(self):
        """builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """
        chain = np.zeros((self.num_particles, 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polymer')
        sig_h = self.ff_reader.get_nbnd_sigma_single_particle('PEG')

        current_pt = np.array([0.0, 0.0, 0.0])
        current_hyb = np.array([0.0*float(poly_bond_length), sig_h/2, 0.0])

        chain[0] = np.array([-1.0, 0.0, 0.0])*poly_bond_length + np.array([1.0, 0.0, 0.0])*poly_bond_length

        for i in range(0, int(self.num_particles/2)):
            ind = i*2 + 1
            pt = np.array([1.0, 0.0, 0.0])*float(poly_bond_length)
            current_pt += pt
            current_hyb += pt
            chain[ind] = current_pt
            chain[ind+1] = current_hyb

        sub_pt = np.array([0, 0, 0])
        for i in range(self.num_particles):
            chain[i] = np.subtract(chain[i], sub_pt)

        return chain

    def get_bond_types(self):
        """utilized to dump_gsd

        :return: list of all the bond types in the system
        """
        return ['polymer'] * self.repeats

    def get_bonds(self):
        """overrides method in LigandAbs to dump GSD files

        :return: n by 2 list of indices that share bond potentials
        """

        bonds = np.zeros((self.repeats, 2))

        bonds[0] = [0, 1]
        for i in range(0, self.repeats-1):
            ind = i*2 + 1
            bonds[i+1] = [ind, ind+2]

        return bonds

    def get_angles(self):
        return []

    def get_dihedrals(self):
        return []
