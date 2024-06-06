"""
:module: UnitsClass
:platform: Unix, Windows
:synopsis: Implements simple class for units

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - updated documentation
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - precisely stated what are construction and simulation units, changed names
..                  - Added documentation
..                  - Introduced angstrom to construction
..                 Alex Travesset <trvsst@ameslab.gov> April 2022
..                  - expanded the range of the class with more unit conversion and change the naming
"""

import numpy as np


class UnitsClass:

    """
    Simple class containing defining units
    """

    def __init__(self, name, lgth, mss, egy, st, bz, pm, dt):
        """The constructor: defines the conversion from construction to simulation units
        
        :param name: units name
        :param lgth: length unit coefficient (conversion from construction to simulation units)
        :param mss:  mass unit coefficient  (conversion from construction to simulation units)
        :param egy: energy unit coefficient (conversion from construction to simulation units)
        :param st: default temp of the simulation in dimensionless units
        :param bz: Boltzmann Constant, in construction units
        :param pm: Vacuum permittivity, in construction units
        :param dt: time step of the simulation, in simulation units
        """

        # this variable is necessary as some nanocrystal files are hard-written in Angstrom.
        self.angstrom_to_construction = None

        # units name
        self.name = name

        # time step in simulation units
        self.dt = dt

        # constants
        self.sim_temp_default = st
        self.boltz = bz
        self.perm = pm

        # units defined
        self.simulation_units_defined = ['length', 'mass', 'energy', 'angle', 'volume', 'moment_inertia', 'time',
                                         'force', 'bond', 'pressure', 'velocity', 'acceleration', 'charge', 'density']

        # basic units
        self.length_construction_to_simulation = lgth
        self.mass_construction_to_simulation = mss
        self.energy_construction_to_simulation = egy

        # from construction to simulation units
        self.angle_construction_to_simulation = 1
        self.volume_construction_to_simulation = lgth**3
        self.moment_inertia_construction_to_simulation = mss*lgth**2
        self.time_construction_to_simulation = lgth*np.sqrt(mss/egy)
        self.force_construction_to_simulation = egy/lgth
        self.bond_construction_to_simulation = egy/lgth**2
        self.pressure_construction_to_simulation = egy/lgth**3
        self.velocity_construction_to_simulation = lgth / self.time_construction_to_simulation
        self.acceleration_construction_to_simulation = lgth / self.time_construction_to_simulation**2
        self.charge_construction_to_simulation = 1/np.sqrt(4*np.pi*self.perm/lgth/egy)
        self.density_construction_to_simulation = mss / self.volume_construction_to_simulation

        # temperature in simulation units
        self.kelvin_temp = self.sim_temp_default / (self.boltz * self.energy_construction_to_simulation)
