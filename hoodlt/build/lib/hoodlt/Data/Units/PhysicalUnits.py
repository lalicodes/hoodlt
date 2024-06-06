"""
:module: PhysicalUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing a physical unit system

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - precisely stated what are construction and simulation units, changed names
..                  - Added documentation
..                 Alex Travesset <trvsst@ameslab.gov> April 2022
..                  - changed the naming of the class variables for easier comprehension
"""

import numpy as np
from hoodlt.Data.Units.UnitsClass import UnitsClass
from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants


class PhysicalUnits(UnitsClass):
    """
    Abstract class which can represent all physical unit systems
    """

    def __init__(self, name, lgth, mss, egy, st, bz, pm, dt, length_to_si, mass_to_si, energy_to_si):
        """

        :param name: name of the units class
        :param lgth: length unit coefficient (conversion from construction to simulation units)
        :param mss: mass unit coefficient  (conversion from construction to simulation units)
        :param egy: energy unit coefficient (conversion from construction to simulation units)
        :param st: default temperature in simulation units
        :param bz: the boltzmann constant in simulation units
        :param pm: the vacuum permittivity constant in simulation units
        :param dt: time step of the simulation, in simulation units
        :param length_to_si: conversion from construction unit of length to SI unit of length (m)
        :param mass_to_si: conversion from construction unit of mass to SI unit of mass (kg)
        :param energy_to_si: conversion from construction unit of energy to SI unit of energy (J)
        """

        super(PhysicalUnits, self).__init__(name, lgth, mss, egy, st, bz, pm, dt)

        # basic conversions to SI
        self.length_construction_to_si = length_to_si
        self.mass_construction_to_si = mass_to_si
        self.energy_construction_to_si = energy_to_si

        # derived conversions to SI
        self.angle_construction_to_si = 1
        self.volume_construction_to_si = length_to_si ** 3
        self.moment_inertia_construction_to_si = mass_to_si * length_to_si ** 2
        self.time_construction_to_si = np.sqrt(mass_to_si * length_to_si ** 2 / energy_to_si)
        self.force_construction_to_si = energy_to_si / length_to_si
        self.bond_construction_to_si = energy_to_si/ length_to_si ** 2
        self.pressure_construction_to_si = energy_to_si / length_to_si ** 3
        self.velocity_construction_to_si = length_to_si / self.time_construction_to_si
        self.acceleration_construction_to_si = length_to_si / self.time_construction_to_si ** 2
        self.charge_construction_to_si = energy_to_si
        self.density_construction_to_si = mass_to_si / self.volume_construction_to_si

        # from simulation to construction units
        str_forward = '_simulation_to_construction'
        str_backward = '_construction_to_simulation'
        for units_name in self.simulation_units_defined:
            obj = getattr(self, units_name + str_backward)
            setattr(self, units_name + str_forward, 1 / obj)

        # from SI to conversion units
        str_forward = '_si_to_construction'
        str_backward = '_construction_to_si'
        for units_name in self.simulation_units_defined:
            obj = getattr(self, units_name + str_backward)
            setattr(self, units_name + str_forward, 1 / obj)

        # from SI to simulation units
        str_forward = '_simulation_to_si'
        str_backward = '_construction_to_si'
        str_next = '_simulation_to_construction'
        for units_name in self.simulation_units_defined:
            obj1 = getattr(self, units_name + str_backward)
            obj2 = getattr(self, units_name + str_next)
            setattr(self, units_name + str_forward, obj1*obj2)

        # from simulation units to si
        str_forward = '_si_to_simulation'
        str_backward = '_simulation_to_si'
        for units_name in self.simulation_units_defined:
            obj = getattr(self, units_name + str_backward)
            setattr(self, units_name + str_forward, 1 / obj)

        # static reference to PhysicalConstants class object
        self.pc = PhysicalConstants

    def relative_temp(self, temp_in_kelvin):
        """
        Gets conversion factor from the units default temperature(self.kelvin_temp) the actual simulation temperature

        :param temp_in_kelvin: temperature of the simulation in Kelvin
        :return: dimensionless value of the temperature
        """

        return temp_in_kelvin/self.kelvin_temp

    def energy_to_kbt(self, temp_in_kelvin):
        """
        Gets conversion factor from values in base energy units to kbT
        :param temp_in_kelvin: temperature of the simulation in Kelvin
        :return: conversion factor from base energy unit to kbT
        """
        return 1 / (self.boltz * temp_in_kelvin)

    def simulation_temp(self, temp_in_kelvin):
        """

        :param temp_in_kelvin: temperature in kelvin
        :return: temperature in simulation units
        """
        return temp_in_kelvin * self.boltz * self.energy_construction_to_simulation

    def sim_frame_to_sec(self, time_step):
        """
        Converts the time of a simulation frame to a time in SI seconds

        :param time_step: the time step in the simulation
        :return: time_step in si units (seconds)
        """

        return time_step * self.dt * self.time_simulation_to_si
