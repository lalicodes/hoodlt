"""
:module: LatticeFunctionalizedConfiguration
:platform: Unix, Windows
:synopsis: create a lattice configuration of NCs

.. moduleauthor:: Xun Zha <xzha@iastate.edu> November 2018
.. history:
..                Xun Zha <xzha@iastate.edu> June 2021
..                  - changed parameter rotation to orientations, removed an unnecessary method
..
..                Alex Travesset <trvsst@ameslab.gov> April 2022
..                  - Enhancements following adoption of hoomd v3
..
..                Alex Travesset <trvsst@ameslab.gov> July 2022
..                  - changed naming so as to allow for lattice constants with decimal numbers
"""

import copy
import numpy as np
import hoodlt.Utils.LatticeNeighbors as Ln
from hoodlt.Data.Modelconfigurations.FunctionalizedConfiguration import FunctionalizedConfiguration


class LatticeFunctionalizedConfiguration(FunctionalizedConfiguration):

    def __init__(self, lat, ncs):
        """Lattice NC Configuration

        :param lat: lattice object
        :param ncs: list of FunctionalizedParticle objects,
        """

        self.lat = lat  # lattice object
        self.neighb = Ln.LatNeighbor(lat)  # lattice neighbor object

        # position vectors for each NC (need to be transposed from [3, num_particles] -> [num_particles, 3])
        vector = np.transpose(self.neighb.l_vec)
        # configuration box [Lx, Ly, Lz, xy, xz, yz]
        a_v = self.neighb.d_d.a_v
        self.box_l = np.concatenate((np.diag(a_v), [a_v[0, 1]/a_v[1, 1], a_v[0, 2]/a_v[2, 2], a_v[1, 2]/a_v[2, 2]]))

        # if only given one NC per type, generate the rest of the NCs for the superlattice
        if len(ncs) == sum(lat.num_pnts()):
            pass
        elif len(ncs) == len(lat.typ):
            ncs = [copy.deepcopy(ncs[self.neighb.typ[i]]) for i in range(len(self.neighb.typ))]
        else:
            raise ValueError('Given number of NCs is inconsistent with lattice!')

        super(LatticeFunctionalizedConfiguration, self).__init__(ncs, vector, self.box_l)

    def add_bonds(self, dict_bonds):
        """

        param dict_bonds: dictionary, key is bond name, value is tuple (particle type, degree=1 (nn), 2 (nnn), etc..)
        """

        # set all lists to zero, forcing all bonds to be updated only once
        self.bonds = []
        self.bonds_types = []
        self.bonds_dist = []

        for itm in dict_bonds.items():
            self.bonds_types.append(itm[0])
            ind_type, degree = itm[1]
            n_list = self.neighb.bonds_in_lattice(ind_type, degree)
            self.bonds.append(n_list)
            num = self.neighb.base[n_list[0][1]]
            self.bonds_dist.append(self.neighb.u_dist[num][degree])

    def get_name(self):
        """gives the name of generated lattice

        :return: name of the generated lattice
        """

        name = ""

        # append nc info
        name += self._build_string(self.particles)

        # append lattice type info
        name += '_c'+self.lat.name()[0].title()

        # append l info
        name += ('_l' + str(self.lat.l[0]))
        if len(self.lat.typ) != 1:
            name += (str(self.lat.l[1]) + str(self.lat.l[2]))

        # append lattice constant
        str_lat = str(self.lat.a_val()).replace('.', 'p')
        name += '_a' + str_lat

        # append solvent information, if there is solvent
        if len(self.solvent) > 0:
            name += '_s' + self._build_string(self.solvent)

        # append substrate info
        if len(self.substrates) > 0:
            name += '_b' + self._build_string(self.substrates)

        # append forcefield info
        name += '_ff' + self.ff_reader.name.title()

        return name
