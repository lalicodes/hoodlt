"""
:module: LigandAbs
:platform: Unix, Windows
:synopsis: Defines the abstract classes used to store ligands

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> April 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - removed a now unecessary function
..                  - Moved many methods to Analysis
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - align moved to BasicSystemEntity
..                  - chain_vector moved to BasicSystemEntity as get_vector
..                  - ligand_name is now set as 'name' in basic entity class
..                Alex Travesset <trvsst@ameslab.gov>  October 2018
..                     - minor extension to accommodate more complex ligands
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                     - added special coulomb and lj functionality
..                Alex Travesset <trvsst@ameslab.gov> May 2022
..                     - simplified class
"""


import numpy as np
from hoodlt.Data.Modelconfigurations.BasicSystemEntity import BasicSystemEntity


class LigandAbs(BasicSystemEntity):
    """
    Defines a general ligand
    """

    def __init__(self, repeats, forcefield, num_particles, name):
        """

        :param repeats: the number of repeat units in the chain
        """

        super(LigandAbs, self).__init__(forcefield, num_particles, name)

        self.repeats = repeats

    def get_vector(self):
        """Computes the length of the chain

        :return: vector from first to last position
        :rtype: numpy pos array [x,y,z]
        """

        return np.subtract(self.position[self.num_particles-1], self.position[0])

    def get_name(self):
        """Creates a name to be used by different files

        :return: string with name
        """
        name = self.name

        if self.repeats > 1:
            name += '-n' + str(self.repeats)

        return name

    def get_bonds(self):
        """
        used by other classes to add bonds when the gsd is dumped
        :return: an n by 2 list of index pairs that are bonded together 
        """

        bonds = np.zeros((self.num_particles-1, 2), dtype=int)
        for i in range(0, self.num_particles-1):
            bonds[i] = np.array([i, i+1], dtype=int)
        return bonds

    def get_angles(self):
        """
        
        used by other classes to add angles when the gsd is dumped
        :return: an n by 3 list of indexes that have an angle together  
        """
        angles = np.zeros((self.num_particles-2, 3), dtype=int)
        for i in range(0, self.num_particles-2):
            angles[i] = np.array([i, i+1, i+2], dtype=int)
        return angles

    def get_dihedrals(self):
        """
        
        used by other classes to add dihedrals when the gsd is dumped
        :return: an n by 4 list of indexes that have a dihedral together 
        """
        dis = np.zeros((self.num_particles-3, 4), dtype=int)
        for i in range(0, self.num_particles-3):
            dis[i] = np.array([int(i), int(i+1), int(i+2), int(i+3)], dtype=int)
        return dis
