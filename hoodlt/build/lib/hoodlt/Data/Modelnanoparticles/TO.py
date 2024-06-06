"""
:module: TO
:platform: Unix, Windows
:synopsis: Defines the class for truncated octahedra

.. moduleauthor:: Xun Zha <xzha@iastate.edu> November 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - rewrote class so it is more general
..                  - reworked class so it fits with BasicSystemEntity
..                  - name is now set in basic entity class
"""

from hoodlt.Data.Modelnanoparticles.CoreWithPositionsInFile import CoreWithPositionsInFile
from hoodlt.Data.Modelnanoparticles.TOs.TO_info import ConfigInfo


class TO(CoreWithPositionsInFile):
    def __init__(self, ff, num_core_particles, version, num_graft_sites, core_types=['Au'], scale=None,
                 name_rigid_center=None):
        """

        :param ff: the name of the forcefield to be used to construct this object
        :param num_core_particles: number of particles on the shell of the core
        :param version: degeneracy
        :param num_graft_sites: number of grafting sites on the core
        :param core_types: type of the core particles list should either be length 1 or length num_core_particles
        :param scale: factor by which to scale all the positions on the core. If left to default, positions and grafting sites will be scaled using the fcc nearest neighbor distances
        :param name_rigid_center: name of the fictional rigid center that will be used for this core, without the typical '_' prefix
        """


        # read configuration information
        info = ConfigInfo(num_core_particles, version, core_types[0])
        self.a_nn = info.a_nn
        self.radius = info.radius
        self.area = info.area
        self.volume = info.volume

        # determine how things will be scaled
        if scale is None:
            if len(list(set(core_types))) > 1:
                raise ValueError("Can't scale values by FCC nearest neighbor distance if there are multiple core types")
            factor = self.a_nn
        else:
            factor = scale

        # make the file name
        name = 'N'+str(int(num_core_particles))+'_v'+str(int(version))+'.txt'
        file_name = "TOs/" + name

        # call super
        super(TO, self).__init__(ff, num_core_particles, file_name, core_types, factor, name_rigid_center)

        # set the grafting sites, since the files don't explicitly have them
        self.graft_num = num_graft_sites
        self.add_grafters()
