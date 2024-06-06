"""
:module: CoreWithGraftingSitesFromFile
:platform: Unix, Windows
:synopsis: Defines an Abstract class which takes both positions and grafting sites from a file in order to construct
a general core

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> July 2019
.. history:
                Alex Travesset <trvsst@ameslab.gov> June 2021
                    - fixed some documentation errors
                    - Reading of the file consistent with the units in the forcefield (adding Angstrom to Construction)
"""

import numpy as np
import importlib_resources
from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs

class CoreWithPositionsInFile(NanoAbs):
    """
    This class is the constructor for any core which takes core particle positions and possibly grafting positions
    from a text file. It should be the abstract class for any core which gets its positions from a file
    """

    def __init__(self, ff, num_core_particles, file_with_positions, core_types=['Au'], scale=1, name_rigid_center=None,
                 num_particles_for_name=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param file_with_positions: file where positions are included
        :param num_core_particles: number of particles on the shell of the core
        :param core_types: type of the core particles list should either be length 1 or length num_core_particles
        :param graft_types: list of types of grafters that will be used to graft to this core
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting
        sites will be in angstroms
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the
        typical '_' prefix
        :param num_particles_for_name: number of particles on the core (includes also grafting sites)
        """

        # check inputs
        types_core = self._check_and_create_core_types(core_types, num_core_particles)

        # get the number of particles right for the name of this core
        if num_particles_for_name is None:
            particles_on_core = num_core_particles
        else:
            particles_on_core = num_particles_for_name

        # create the name for the core's rigid center
        dist_core_types = list(set(types_core))
        if name_rigid_center is None:
            name = self._make_string(dist_core_types) + str(particles_on_core)
        else:
            name = name_rigid_center + str(particles_on_core)

        # once you have the name, call the super constructor
        super(CoreWithPositionsInFile, self).__init__(ff, num_core_particles+1, name)

        # set the masses and charges
        for i, typ in enumerate(types_core):
            self.mass[1+i] = self.ff_reader.get_molecular_weight(typ)
            self.charge[1+i] = self.ff_reader.get_charge(typ)
        self.mass[0] = np.sum(self.mass[1:])

        # type data
        self.types = ['_'+name] + dist_core_types
        self.typeid = [self.types[0]] + types_core

        #get the units of the forcefield
        angstrom_to_construction = self.ff_reader.get_units().angstrom_to_construction

        # particle positions and grafting sites
        d_name = 'Data/Modelnanoparticles/' + file_with_positions
        ref = importlib_resources.files('hoodlt') / d_name
        with importlib_resources.as_file(ref) as path:
            positions = (np.genfromtxt(path) * scale).tolist()
        self.position = np.array([[0.0, 0.0, 0.0]] + positions[:num_core_particles])*angstrom_to_construction

        # grafting info
        self.graft_sites = np.array(positions[num_core_particles:])*angstrom_to_construction
        self.graft_num = len(self.graft_sites)

        # moment of inertia in rest frame
        self.moment_inertia[0] = np.diag(self.moment_of_inertia())

    @staticmethod
    def _make_string(list):
        """
        Makes a string out the the elements in the list

        :param list: the input list
        :return:
        """

        string = ""
        for elt in list:
            string += str(elt) + "-"

        return string[:-1]

    @staticmethod
    def _check_and_create_core_types(list_core_types, num_core_particles):
        """
        Checks to make sure the list of core types is valid based on the number of core particles, and returns an
        appropriate list of core types for later use.

        :param list_core_types: list of types on the core
        :param num_core_particles: number of particles on the core
        :return: a list of core types for use later in the constructor
        """

        if len(list_core_types) == 1:
            types_core = [list_core_types[0] for _ in range(num_core_particles)]
        elif len(list_core_types) != num_core_particles:
            raise ValueError("Length of the core_types list should be either 1 or " + str(num_core_particles))
        else:
            types_core = list_core_types

        return types_core
