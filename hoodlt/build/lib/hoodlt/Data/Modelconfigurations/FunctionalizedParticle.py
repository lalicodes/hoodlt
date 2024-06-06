"""
:module: FunctionalizedParticle
:platform: Unix, Windows
:synopsis: combines a NanoAbs implementation with a LigandAbs implementation 

.. moduleauthor:: Curt Waltmann <waltmann@iastate.edu> April 2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - removed an unecessary function
..                  - updated the grafting scheme, cores can now dump all their positions.
..                  - Moved many methods that should be in analysis to Analysis.Analyze
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - calls SystemEntity constructor to set the forcefield reader
..                  - altered max_radius and average_hydrodynamic_radius to account for ncs with images
..                  - altered remove and apply units slightly
..                  - made to inherit from SystemEntity class
..                  - added dumping ncs with impropers
..                  - updated code for orienting/de-orienting ncs
..                  - removed some unecessary methods
..                  - removed references to body_count variable when dumping
..                  - removed the '_' in front of a few method names
..                  - Deprecated the chain maps
..                Tommy Waltmann <tomwalt@iastate.edu> May 2019
..                  - removed unecessary fields stored in the constructor
..                  - made sure the ligands are grafted to the core surface as soon as the object is constructed
..                  - rewrote many member functions to be simpler, easier to understand
..                Alex Travesset <trvsst@ameslab.gov>  October 2018
..                     - made class closed for modification
..                Xun Zha <xzha@iastate.edu>  January 2018
..                     - changed condition `if i != self.core.graft_type' to `if i not in self.core.graft_type'
..                     - added condition `if exist attribute mass_total, mass_center=mass_total'
..                Jacob Austin <jaustin2@iastate.edu> April 2020
                    - added special coulomb and lj functionality
..                Xun Zha <xzha@iastate.edu> Sept 2020
..                  - changed `nano_core.graft_sum' to `np.sum(nano_core.graft_num)'
..                Alex Travesset <trvsst@ameslab.gov> June 2022
..                  - eliminated the dependence on hoomd classes
..                  - Made it consistent with HOOMD v3
"""

import copy as cp
import numpy as np
from hoodlt.Data.Modelconfigurations.SystemEntity import SystemEntity
from hoodlt.Analysis.Analyze.AnalyzeNanoparticle import AnalyzeNanoparticle as AnalyzeNC
from hoodlt.Data.Modelconfigurations.FunctionalizedConfiguration import ParticleToConfiguration
from hoodlt.Groups.Quat import Quat


class FunctionalizedParticle(SystemEntity):
    """
    Defines a nanoparticle functionalized with ligands. For documentation of methods that don't have any, look in the
    SystemEntity class
    """

    def __init__(self, nano_core, ligand):
        """

        :param nano_core: NanoAbs object
        :param ligand: LigandAbs array
        """

        # make sure the number of ligands is correct
        if np.sum(nano_core.graft_num) != len(ligand):
            raise ValueError("You are trying to graft " + str(len(ligand)) + " ligands to a core which has " +
                str(nano_core.graft_num) + " grafting sites")

        self.core = nano_core
        self.ligands = ligand

        # assumes core and ligands read from the same forcefield
        super(FunctionalizedParticle, self).__init__(cp.deepcopy(self.core.ff_reader))

        # graft the ligands to the surface of the NC
        for i, pos in enumerate(self.core.graft_sites):
            self.ligands[i].align(self.ligands[i].position[0], pos - self.core.position[0])
            self.ligands[i].shift(-1 * self.ligands[i].position[0])
            self.ligands[i].shift(pos)

        # update the moment of inertia for the core, now that the grafters are part of the rigid body
        self._update_moment_inertia()

    def _update_moment_inertia(self):
        """
        Updates the moment of inertia and total mass for the core, once the ligands are grafted to its surface

        :return: None
        """

        # add the masses and positions to a copy of the original core object
        # this is sloppy coding, but it gets the job done
        core_copy = cp.deepcopy(self.core)

        mass_grafters = [lig.mass[0] for lig in self.ligands]
        core_copy.mass = np.array(list(self.core.mass) + mass_grafters)
        core_copy.position = np.array(list(self.core.position) + list(self.core.graft_sites))

        self.core.moment_inertia[0] = np.diag(core_copy.moment_of_inertia())
        self.core.mass[0] += np.sum(mass_grafters)

    def get_num_particles(self):

        total = self.core.num_particles
        for lig in self.ligands:
            total += lig.num_particles

        return total

    def add_to_snap(self, snap, tag):

        nc_body_index = tag

        tag = self.core.add_to_snap(snap, tag)

        for lig in self.ligands:
            tag_before = tag
            tag = lig.add_to_snap(snap, tag)

            # grafter on the ligand must be part of the core's rigid body
            snap.particles.body[tag_before] = nc_body_index

        return tag

    def get_distinct_bonds_data(self, attr):

        a_types = getattr(self.core, 'get_distinct_bonds_data')(attr)
        for lig in self.ligands:
            for a_typ in getattr(lig, 'get_distinct_bonds_data')(attr):
                if a_typ not in a_types:
                    a_types.append(a_typ)
        return a_types

    def get_distinct_particle_types(self):

        p_types = []
        for i in self.core.types:
            p_types.append(i)
        for lig in self.ligands:
            for x in lig.types:
                if x not in p_types:
                    p_types.append(x)
        return p_types

    def _get_box_for_dumping(self):

        mat = np.array([1000, 1000, 1000, 0, 0, 0])
        ac = AnalyzeNC(self, mat)
        l_val = int(2*ac.max_radius() + 25)
        return np.array([l_val, l_val, l_val, 0.0, 0.0, 0.0])

    def shift(self, vector):
        """
        shifts the entire particle by some vector
        :param vector: vector the particle is shifted by [x,y,z]
        :return: nothing
        """

        self.core.shift(vector)
        for lig in self.ligands:
            lig.shift(vector)

    def apply_units(self):
        """
        Apply units
        """

        self.core.apply_units()

        for lig in self.ligands:
            lig.apply_units()

    def rotate(self, quat):
        """
        rotate the functionalized particle

        :param quat: quternion (as a 4-tuple)
        """
        for lig in self.ligands:
            lig.rotate(self.core.position[0], quat)
        self.core.rotate(quat)

    def rotate_actual(self, quat):
        """
        rotate the functionalized particle (this is used in drawing rotated nanoparticles)

        :param quat: quternion (as a 4-tuple)
        """

        for lig in self.ligands:
            lig.rotate(self.core.position[0], quat)
        self.core.rotate_actual(quat)

    def body_frame(self):
        """
        returns the coordinates of the body frame
        """
        qt = Quat()
        quat = tuple(qt.inverse(self.core.orientation[0]))
        self.rotate(quat)
        self.shift(-self.core.position[0])

    def scale(self, val):
        """
        takes all the positions and multiples them by the val; essentially rescaling the system
        :param val: scaling factor
        :return: void
        """

        self.core.scale(val)
        for lig in self.ligands:
            lig.scale(val)

    def get_name(self):

        ligs = self._get_distinct_ligands()

        # build string from all the distinct ligands
        lig_string = ""
        for lig in ligs:
            lig_string += lig.get_name() + "-"

        return self.core.get_name() + '-' + lig_string[:-1]

    def _get_distinct_ligands(self):
        """
        Returns a list of distinct ligands on this core

        :return: a list of LigandAbs objects
        """

        list_distinct_ligs = []
        for lig in self.ligands:
            # determine if lig is distinct from the ligs in list_distinct_ligs
            distinct = True
            for dist_lig in list_distinct_ligs:
                if lig == dist_lig:
                    distinct = False
                    break
            if distinct:
                list_distinct_ligs.append(lig)

        return list_distinct_ligs

    def configuration(self):
        """

        :return:the particle as a FunctionalizedConfiguration
        """

        return ParticleToConfiguration(self)

    def __eq__(self, other):
        """
        Compares two objects to determine if they are equal. We define one nanoparticle to be equal to another if they
        have the same core, and ligand objects. Assumes 1 ligand object per nanoparticle
        :param other: other FunctionalizedParticle object
        :return:
        """
        return self.core.name == other.core.name and \
               self.ligands[0].name == other.ligands[0].name and \
               self.ligands[0].repeats == other.ligands[0].repeats
