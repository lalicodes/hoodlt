"""
:module: GsdParser
:platform: Unix
:synopsis: Class which parses gsd files to acquire information for setting up and running simulations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - rigid body relative positions are now set in the rest frame of the rigid bodies
..                  - Added documentation
..
..                Alex Travesset <trvsst@ameslab.gov>, April 2022
..                  - Added additional functions to identify nanoparticle tags
"""


from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions


class GsdParser:
    """
    Class which parses gsd files for information needed to set up hoomd simulations
    """

    def __init__(self, snap):
        """

        :param snap: hoomd snapshot object which represents the initial state of a simulation
        """

        self.snap = snap
        self.ctr_bond_prefix = 'CTR-CTR'
        l_b = self.snap.configuration.box
        self.box = BoxFunctions(l_b)

    def has_charges(self):
        """
        Determines whether or not the system has non-zero charges

        :return: True if the system has non-zero charges, false otherwise
        """

        for i in range(self.snap.particles.N):
            if self.snap.particles.charge[i] != 0.0:
                return True

        return False

    def has_rigid_bodies(self):
        """
        Determines whether or not the system has rigid bodies

        :return: True if the system has rigid bodies, false otherwise
        """

        for i in range(self.snap.particles.N):
            if self.snap.particles.body[i] != -1:
                return True

        return False

    def integration_group(self):
        """
        Returns an index in a list defined in HoomdSimulation, tells the simulation which group to integrate over.
        1 -> integrate over union of rigid centers and non-rigid particles
        0 -> integrate over all particles

        :return: an index in a list defined in HoomdSimulation constructor
        """

        if self.has_rigid_bodies():
            return 1
        else:
            return 0
