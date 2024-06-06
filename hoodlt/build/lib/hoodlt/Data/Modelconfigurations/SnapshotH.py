"""
:module: SnapshotH
:platform: Unix, Windows
:synopsis: class that defines the properties of a HOOMD snapshot
.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov> March 2022
.. history:
..
..                Alex Travesset <trvsst@ameslab.gov> March 2022
                    - defined the class to be used with hoomd v3.0 and above
"""

import copy as cp


class SnapshotDefault(object):

    def __init__(self):
        """
        creates an empty snapshot with default to zero
        """

        self.N = 'N'
        self.types = 'types'
        self.box = 'box'

        # set all quantities to be stored to their defaults, subclasses will override some of these quantities
        self.particles_default = {'position': [0.0, 0.0, 0.0], 'image': [0, 0, 0], 'typeid': 0,
                                  'velocity': [0.0, 0.0, 0.0], 'acceleration': [0.0, 0.0, 0.0], 'charge': 0.0,
                                  'mass': 1.0, 'diameter': 1.0, 'body': -1, 'orientation': [1.0, 0.0, 0.0, 0.0],
                                  'moment_inertia': [0.0, 0.0, 0.0]}

        self.particles_units = {'position': 'length', 'image': None, 'typeid': None, 'velocity': 'velocity',
                                'acceleration': 'acceleration', 'charge': 'charge', 'mass': 'mass',
                                'diameter': 'length', 'body': None, 'orientation': None,
                                'moment_inertia': 'moment_inertia'}

        self.bonds_data = ['bonds', 'angles', 'dihedrals', 'impropers', 'pairs']
        dict_empty = {'types': [], 'typeid': [], 'group': []}
        for bd in self.bonds_data:
            setattr(self, bd+'_default', cp.deepcopy(dict_empty))

        self.constraints_data = ['constraints']
        dict_empty = {'value': [], 'group': []}
        for ct in self.constraints_data:
            setattr(self, ct + '_default', cp.deepcopy(dict_empty))
