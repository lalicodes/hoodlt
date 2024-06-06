"""
:module: AuS
:platform: Unix, Windows
:synopsis: Defines the class for gold truncated octahedra grafted with sulfur ligand head

.. moduleauthor:: Xun Zha <xzha@iastate.edu> Dec 2020
.. history:
                  Xun Zha <xzha@iastate.edu> June 2021
                    - Reading of the file consistent with the units in the forcefield (adding Angstrom to Construction)

                  Alex Travesset <trvsst@ameslab.gov> October 2023
                    - Used importlib.resources as pkk_resources is outdated
"""

import numpy as np
import importlib_resources
from hoodlt.Data.Modelnanoparticles.NanoAbs import NanoAbs


class AuS(NanoAbs):
    def __init__(self, forcefield, num_pnts, version):
        """

        :param forcefield: the name of the forcefield to be used to construct this object
        :param num_pnts: number of Au particles in the core
        :param version: degeneracy of num_pnts
        """

        # check the range for num_pnts and version
        if num_pnts < 38 or num_pnts > 4033:
            raise ValueError('Parameter num_pnts needs to be between 38 and 4033')
        if version not in [1, 2, 3]:
            raise ValueError('Parameter version needs to be 1 or 2 or 3')

        # read configuration information
        d_name = 'Data/Modelnanoparticles/AuSs/'+str(int(num_pnts)) + '_' + str(int(version)) + '_AuS.txt'
        ref = importlib_resources.files('hoodlt')/d_name
        with importlib_resources.as_file(ref) as path:
            data = np.loadtxt(path).T
        typeid = data[0]
        pos = data[1:].T

        # number of gold atoms of the bulk and number of gold atoms on the surface
        num_bulk_particles = int(num_pnts)
        num_particles = sum(typeid == 0)

        # name of the core rigid center, call the super constructor
        name = 'cAu'+str(num_bulk_particles)+'S'
        super(AuS, self).__init__(forcefield, 1+num_particles, name)

        # get the units of the forcefield
        angstrom_to_construction = self.ff_reader.get_units().angstrom_to_construction

        # particle positions and grafting sites
        self.position = np.vstack((np.zeros([1, 3]), pos[typeid == 0])) * angstrom_to_construction
        self.graft_sites = pos[typeid == 1] * angstrom_to_construction
        self.graft_num = len(self.graft_sites)

        # type data
        self.types = ['_'+name] + ['Au']
        self.typeid = [self.types[0]] + [self.types[1]] * num_particles

        # masses and charges
        m_au = self.ff_reader.get_molecular_weight('Au')
        c_au = self.ff_reader.get_charge('Au')
        for i in range(num_particles):
            self.mass[1+i] = m_au
            self.charge[1+i] = c_au
        self.mass[0] = m_au * num_particles

        # moment of inertia in rest frame
        self.moment_inertia[0] = np.diag(self.moment_of_inertia())
