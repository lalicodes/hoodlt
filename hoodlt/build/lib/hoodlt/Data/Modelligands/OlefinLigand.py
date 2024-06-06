"""
:module: OlefinLigand
:platform: Unix, Windows
:synopsis: Implements the olefin version of the abstract ligand class 

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> September 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
..
..                Alex Travesset <trvsst@ameslab.gov> June 2021
..                  - Fixed previous issues with the constructor
..                  - Documentation issues
"""

from __future__ import division
import numpy as np
from hoodlt.Data.Modelligands.LigandAbs import LigandAbs


class OlefinLigand(LigandAbs):
    """
    Defines the Olefin ligand
    """

    def __init__(self, carbons, ff, double_locations, cis):
        """

        :param carbons: number of carbons in the chain. Total length will have one more, S
        :param double_locations: location of the first CH in the double bond, must be at least 2 carbons apart
        :param cis: list of booleans length of double_locations; True is cis; Trans is False 
        """

        print('This class is not working at the moment')

        for ind, loc in enumerate(double_locations):
            for ind2 in range(ind + 1, len(double_locations)):
                if np.abs(loc - double_locations[ind2]) < 3:
                    raise ValueError('All double bonds must be separated by at least 2 carbons')

        if not len(cis) == len(double_locations):
            raise ValueError('must have as many values for cis as double bonds')

        # ligand name
        ligand_name = 'Olefin'
        for ind, spot in enumerate(double_locations):
            ligand_name += '-'
            ligand_name += str(spot)
            ligand_name += '-'
            if cis[ind]:
                ligand_name += 'Cis'
            else:
                ligand_name += 'Trans'

        # 3rd argument may be wrong, I don't know how this class is supposed to work
        super(OlefinLigand, self).__init__(carbons, ff, len(double_locations)*2 + carbons - 1, ligand_name)

        # cis = [x for _, x in sorted(zip(double_locations, cis))]
        # double_locations.sort(key=float, reverse=True)
        # number of particles in the chain
        self.typeid = ['S']
        last = -1
        for i in double_locations:
            self.typeid += ['CH2'] * (i-last-2)
            self.typeid += ['CH'] * 2
            last = i
        self.typeid += ['CH2'] * (carbons-last-2)
        self.typeid += ['CH3']

        weights = [self.ff_reader.get_molecular_weight('S'),
                   self.ff_reader.get_molecular_weight('CH'),
                   self.ff_reader.get_molecular_weight('CH2'),
                   self.ff_reader.get_molecular_weight('CH3')]

        self.types = ['S', 'CH', 'CH2', 'CH3']

        self.mass = [weights[self.types.index(typ)] for typ in self.typeid]

        # particle positions (Angstroms)
        self.position = self.build_chain(cis)

        # ligand name
        self.ligand_name = 'Olefin'
        for ind,spot in enumerate(double_locations):
            self.ligand_name += '-'
            self.ligand_name += str(spot)
            self.ligand_name += '-'
            if cis[ind]:
                self.ligand_name += 'Cis'
            else:
                self.ligand_name += 'Trans'

        self.cis = cis

    def build_chain(self, cis):
        """
        called by initializer to get the positions

        :return: the position array for the hydrocarbon chain
        """

        chain = np.zeros((np.sum(self.num_particles), 3))

        single_bond_length = self.ff_reader.get_bond_r0('CH2-CH2')
        double_bond_length = self.ff_reader.get_bond_r0('CH-CH')

        double_angle = np.deg2rad(119)

        s_bond_length = self.ff_reader.get_bond_r0('S-CH2')
        c_angle = self.ff_reader.get_angle_t0('CH2-CH2-CH2')
        s_angle = self.ff_reader.get_angle_t0('S-CH2-CH2')

        ch2vec1 = np.array([single_bond_length * np.sin(c_angle / 2.0), single_bond_length * np.cos(c_angle / 2.0), 0.0])
        ch2vec2 = np.array([single_bond_length * np.sin(c_angle / 2.0), -single_bond_length * np.cos(c_angle / 2.0), 0.0])
        ch2vec3 = np.array([-single_bond_length * np.sin(c_angle / 2.0), single_bond_length * np.cos(c_angle / 2.0), 0.0])
        ch2vec4 = np.array([-single_bond_length * np.sin(c_angle / 2.0), -single_bond_length * np.cos(c_angle / 2.0), 0.0])

        double_vec_1 = np.array([double_bond_length * np.sin(double_angle / 2.0), 0,double_bond_length * np.cos(double_angle / 2.0)])
        double_vec_2 = np.array([double_bond_length * np.sin(double_angle / 2.0), 0,- double_bond_length * np.cos(double_angle / 2.0)])

        s_pos = np.array([-s_bond_length * np.sin(s_angle / 2.0), s_bond_length * np.cos(s_angle / 2.0), 0.0])
        #print(c_angle)
        #print(double_angle)

        current = np.array([0.0, 0.0, 0.0])
        chain[0] = s_pos
        chain[1] = current

        count = 0
        self.typeid.append('CH4')

        for i in range(2, self.repeats + 1):
            if i % 2 == 0:
                if self.typeid[i] == 'CH2' or self.typeid[i] == 'CH3':
                    if self.typeid[i-1] == 'CH' and self.typeid[i-2] == 'CH':
                        if cis[count]:
                            current = np.add(current, ch2vec4)
                        else:
                            current = np.add(current, ch2vec1)
                        count += 1
                    else:
                        current = np.add(current, ch2vec1)
                elif self.typeid[i] == 'CH' and self.typeid[i-1] == 'CH':
                    current = np.add(current, double_vec_1)
                else:
                    current = np.add(current, ch2vec1)
            else:
                if self.typeid[i] == 'CH2' or self.typeid[i] == 'CH3':
                    if self.typeid[i-1] == 'CH' and self.typeid[i-2] == 'CH':
                        if cis[count]:
                            current = np.add(current, ch2vec3)
                        else:
                            current = np.add(current, ch2vec2)
                        count += 1
                    else:
                        current = np.add(current, ch2vec2)
                elif self.typeid[i] == 'CH' and self.typeid[i-1] == 'CH':
                    current = np.add(current, double_vec_2)
                else:
                    current = np.add(current, ch2vec2)

            chain[i] = current
        self.typeid.remove('CH4')

        for i in range(self.repeats + 1):
            chain[i] = np.subtract(chain[i], s_pos)

        return chain

    def get_angle_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd

        :return: list of all the angle types in the chain 
        """

        angles = ['S-CH2-CH2']
        for i in range(1, self.repeats-1):
            if 'CH' in self.typeid[i:i+3]:
                angles.append('CH-CH2-CH2')
            else:
                angles.append('CH2-CH2-CH2')
        return angles

    def get_bond_types(self):
        """
        override of the method in LigandAbs
        used to dump_gsd
        :return: list of all the bond types in the system 
        """
        bonds = ['S-CH2']
        for i in range(1, self.repeats):
            if 'CH' == list(set(self.typeid[i:i+2]))[0]:
                bonds.append('CH-CH')
            else:
                bonds.append('CH2-CH2')
        return bonds

    def get_dihedral_types(self):
        """
        :return: list of all the dihedral types in the system
        """
        dihedrals = ['JCH2-CH2-CH2-CH2']

        count = 0

        for i in range(1, self.repeats - 2):
            if 'CH' == self.typeid[i:i+4][1] and 'CH' == self.typeid[i:i+4][2]:
                #print(i)
                #print(self.typeid[i:i+4])
                if self.cis[count]:
                    dihedrals.append('CIS')
                else:
                    dihedrals.append('TRANS')
                count += 1
            elif 'CH' == self.typeid[i:i+4][0] and 'CH' == self.typeid[i:i+4][1] or 'CH' == self.typeid[i:i+4][2] and 'CH' == self.typeid[i:i+4][3] :
                dihedrals.append('1-BUTENE')
            else:
                dihedrals.append('CH2-CH2-CH2-CH2')

        return dihedrals
