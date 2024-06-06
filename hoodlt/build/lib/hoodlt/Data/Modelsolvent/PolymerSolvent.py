"""
:module: PolymerSolvent
:platform: Unix, Windows
:synopsis: Implements the generic polymer solvent version of the abstract ligand class

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu> May 2018
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

class PolymerSolvent(SolventAbs):
    """Defines a solvent consisting of a linear polymer with hybridizable flanking sites
    """

    def __init__(self, repeats, ff):
        """

        :param repeats: the number of repeats on the solvent
        :param ff: the name of the forcefield used to create this solvent
        """

        super(PolymerSolvent, self).__init__(repeats, ff, repeats, 'PolymerSolvent')

        for i in range(repeats):
            self.mass[i] = self.ff_reader.get_molecular_weight('O')

        # types
        self.types = ['O']
        self.typeid = ['O']*self.repeats

        # positions
        self.position = self.build_chain()

    def build_chain(self):
        """
        builds polymer chain positions based on number in chain

        :return: the position array for the chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        poly_bond_length = self.ff_reader.get_bond_r0('polybond', qualifier='fene')

        current_pt = np.array([0.0, 0.0, 0.0])

        for i in range(self.repeats):

            pt = np.array([1.0, 0.0, 0.0])*float(poly_bond_length)*2/3

            current_pt += pt

            chain[i] = current_pt

        sub_pt = np.array([0, 0, 0])
        for i in range(self.num_particles):
            chain[i] = np.subtract(chain[i], sub_pt)

        return chain

    def get_bond_types(self):
        """
        utilized to dump_gsd

        :return: list of all the bond types in the system
        """

        return ['polybond'] * (self.repeats-1)

    def get_bonds(self):
        """
        overrides method in BasicSystemEntity to dump GSD files

        :return: n by 3 list of indices that share angle potentials
        """

        bonds = np.zeros((self.num_particles-1, 2))

        for i in range(self.repeats-1):
            bonds[i] = [i, i+1]

        return bonds

    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[-1] - self.position[0]
