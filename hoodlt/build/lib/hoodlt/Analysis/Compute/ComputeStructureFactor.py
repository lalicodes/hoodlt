"""
:module: ComputeStructureFactor
:platform: Unix, Windows
:synopsis: computes the structure factor of a lattice system

.. moduleauthor: Xun Zha <xzha@iastate.edu>, November 2018
.. history:
..               Alex Travesset <trvsst@ameslab.gov>, August 2022
..                  - rewrote the class to make it compatible with all new changes
"""


import numpy as np
import numpy.linalg as la
import gsd.hoomd

from hoodlt.Analysis.Collect.CollectHistData import CollectHistData
from hoodlt.Data.Modelconfigurations.Saver import load_config
import hoodlt.Data.Modelconfigurations.MapSnapshot as Ma
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions
import hoodlt.Lattices.reciprocal_vectors as rv


class ComputeStructureFactor(CollectHistData):

    def __init__(self, file_names, num_frames=None):
        """
        initializer for the class

        :param file_names: list of filenames
        :param num_frames: number of frames to include
        """

        # get the write files
        self.file_traj = [gsd.hoomd.open(fn + '_write.gsd', 'rb') for fn in file_names]
        end_val = [len(fl) for fl in self.file_traj]

        super(ComputeStructureFactor, self).__init__(file_names, num_frames=num_frames, end_val=end_val)

        # lattice object and number of (nano)particles in the lattice
        conf = load_config(self.file_names[0])
        self.num_lattice_points = conf.lat.num_pnts()
        snap = self.file_traj[0][0]
        map_snap = Ma.MapSnapshot(conf, snap).list_of_entities
        if len(map_snap) != np.sum(self.num_lattice_points):
            raise ValueError('there is a mistmatch between number of entities and lattice object')
        center_tags = np.zeros([np.sum(self.num_lattice_points)], dtype=int)
        self.lattice_np_type = np.sum(self.num_lattice_points)*[0]
        for ind, dict_p in enumerate(map_snap):
            center_tags[ind] = dict_p['tag-center']
            self.lattice_np_type[ind] = dict_p['tag-center-name']
        self.lattice_ideal_pos = np.zeros([self.num_files, np.sum(self.num_lattice_points), 3])
        self.box_dim = self.num_files*[0]
        for ind, names in enumerate(self.file_names):
            conf = load_config(names)
            self.box_dim[ind] = BoxFunctions(conf.box)
            self.lat = conf.lat
            # volume per particle
            self.lattice_ideal_pos[ind] = np.transpose(conf.neighb.l_vec)
        self.lattice_actual_pos = np.zeros([self.num_files, self.num_frames, np.sum(self.num_lattice_points), 3])

        # get the center positions
        for ind_f, fl in enumerate(self.file_traj):
            ini = self.ini_val[ind_f]
            end = self.end_val[ind_f]
            for ind_t, snap in enumerate(fl[ini:end]):
                pos = snap.particles.position[center_tags]
                self.lattice_actual_pos[ind_f, ind_t] = pos

        # compute vibrations
        self.vbr = np.zeros_like(self.lattice_actual_pos)
        for ind1 in range(self.num_files):
            for ind2 in range(np.sum(self.num_lattice_points)):
                self.vbr[ind1, :, ind2] = self.box_dim[ind1].wrap(self.lattice_actual_pos[ind1, :, ind2, :]
                                                                  - self.lattice_ideal_pos[ind1, np.newaxis, ind2, :])

        self.d_lattice_types = list(set(self.lattice_np_type))
        self.lattice_types = [self.d_lattice_types.index(val) for val in self.lattice_np_type]

    def compute_bfactor(self):
        """
        Computes B-factors

        :return : B-factors for each lattice and particle
        """

        c_fac = 8*np.pi**2/3

        return c_fac*np.average(la.norm(self.vbr, axis=-1)**2, axis=1)

    def compute_avg_displacements(self):
        """
        Computes the average quantity :math:`\langle {\vec u} \rangle`

        :return : average for each lattice, particle and x,y,z
        """

        return np.average(self.vbr, axis=1)

    def compute_lindemann(self):
        """"
        Computes Lindemann criteria

        :return: lindemann criteria for each lattice
        """

        dist = np.amin(self.bond_r0, axis=1)
        c_fac = 3/(8 * np.pi ** 2)
        u_2 = c_fac*self.compute_bfactor()

        return np.sqrt(u_2)/dist

    def compute_avg_lindemann(self):
        """"
        Computes average Lindemann

        :return: average lindemann criteria with error
        """

        mat = self.compute_lindemann()
        res = np.zeros([mat.shape[0], 2])
        res[:, 0] = np.average(mat, axis=1)
        res[:, 1] = np.std(mat, axis=1)

        return res
