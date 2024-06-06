"""
:module: PolymerSolventInternalHyb
:platform: Unix, Windows
:synopsis: Implements the polymer version of the abstract ligand class for use in BINARY hydrogen-bonding simulations

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> October 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Removed unecessary functions
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - edited get_vector() method
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class PolymerSolventInternalHyb(SolventAbs):
    """
    Defines the Polymer Ligand with hybridizable flanking sites
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: the number of repeats on the chain
        :param ff: the name of the forcefield to be used to construct this solvent
        """

        super(PolymerSolventInternalHyb, self).__init__(repeats, ff, 2*repeats, 'PolymerSolventInternalHyb')

        # masses, rigid stuff
        for i in range(repeats):
            ind = i*2
            self.mass[ind] = self.ff_reader.get_molecular_weight('cPAA')
            self.mass[ind + 1] = self.ff_reader.get_molecular_weight('PAA')
            self.moment_inertia[ind] = np.array([0, 0.625, 0.625])
            self.body[ind:ind+2] = ind

        # types
        self.types = ['cPAA', 'PAA']
        self.typeid = ['cPAA', 'PAA']*repeats

        # positions
        self.position = self.build_chain()
        self.shift(-1*self.center_of_mass())

    def build_chain(self):
        """
        builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """

        chain = np.zeros((self.num_particles, 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polymer')
        sig_h = self.ff_reader.get_nbnd_sigma_single_particle('PAA')

        current_pt = np.array([0.0, 0.0, 0.0])
        current_hyb = np.array([0.0*float(poly_bond_length), sig_h/2, 0.0])

        chain[0] = np.array([-1.0, 0.0, 0.0])*poly_bond_length + np.array([1.0, 0.0, 0.0])*poly_bond_length

        for i in range(0, int(self.num_particles/2)):
            ind = i*2
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
        """
        utilized to dump_gsd

        :return: list of all the bond types in the system
        """

        return ['polymer'] * (self.repeats - 1)

    def get_bonds(self):
        """
        overrides method in SolventAbs to dump GSD files

        :return: n by 2 list of indices that share bond potentials
        """

        bonds = np.zeros((self.repeats - 1, 2))

        for i in range(0, self.repeats - 1):
            ind = i*2
            bonds[i] = [ind, ind+2]

        return bonds

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[-2] - self.position[0]
