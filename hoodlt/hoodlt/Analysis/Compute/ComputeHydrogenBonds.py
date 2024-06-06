"""
:module: ComputeHydrogenBonds
:Platform: Windows, Unix
:synopsis: Computes Hydrogen Bonds for a given configuration

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu>, August 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, June 2022
..                  - changed the algorithm
..                  - made it compatible with hoomd v3
<<<<<<< HEAD
..
..                Elizabeth Macias <emacias@iastate.edu>, October 2022
..                  - added hydrogen bonds per peo analysis
=======
>>>>>>> 9b9b57136be69938520d07ecd2db26590be1dd2e
"""
import copy

import numpy as np
import gsd.hoomd
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
import hoodlt.Data.Modelconfigurations.Saver as Sv


class ComputeHydrogenBonds:
    """
    compute hydrogen bonding
    """
    def __init__(self, file_name, num_frames=None, donor_list=None, proton_list=None, acceptor_list=None, eps=1e-3):
        """

        Hydrogen bond calculation following the criteria of Berndt et. al. (1993)

        :param file_name: file to analyze the trajectory
        :param num_frames: number of frames to consider
        :param donor_list: list of donors
        :param proton_list: list of protons
        :param acceptor_list: list of acceptors
        :param eps: precision of the donor-proton constrain distance
        """

        # initialize parameters
        if donor_list is None:
            donor_list = ['OW']
        if proton_list is None:
            proton_list = ['H']
        if acceptor_list is None:
            acceptor_list = ['OW']

        self.donor_list = donor_list
        self.proton_list = proton_list
        self.acceptor_list = acceptor_list

        self.num_hb_types = len(self.donor_list)

        if len(self.proton_list) != self.num_hb_types or len(self.acceptor_list) != self.num_hb_types:
            raise ValueError('provided lists of acceptors, protons and donors are inconsistent')

        self.eps = eps

        # read files
        self.file_name = file_name + '_write.gsd'
        self.sys = gsd.hoomd.open(self.file_name, mode='r')
        conf = Sv.load_config(file_name)
        self.ff = ForceFieldReader(conf.ff_reader.name)
        self.units = self.ff.get_units()

        # frames to compute hbs
        if num_frames is None:
            num_frames = len(self.sys)
        if num_frames > len(self.sys) or num_frames < 0:
            raise ValueError("frame out of bounds")
                                                                     
        self.end_frame = len(self.sys)
        self.ini_frame = self.end_frame-num_frames
        self.num_frames = num_frames

        # hydrogen bond definitions
        r_cut_angstrom = 2.5
        self.r_cut = r_cut_angstrom*self.units.angstrom_to_construction*self.units.length_construction_to_simulation
        self.theta_cut = np.deg2rad(150)
        self.min_cos = -1.0
        self.max_cos = np.cos(self.theta_cut)

        # store hb information
        self.hbonds = self.num_hb_types*[self.num_frames*[0]]
        self.num_adp = self.num_hb_types*[0]

        for ind in range(self.num_hb_types):
            hb_tags, all_tags = self._get_hbonds(self.donor_list[ind], self.proton_list[ind], self.acceptor_list[ind])
            self.hbonds[ind] = hb_tags
            self.num_adp[ind] = all_tags

    def average_mol_hb(self):
        """
        average number of hb per acceptor, donor for each frame

        :return: a numpy array
        """

        return np.sum(np.arange(6)*self.average_hist_hb(), axis=2)

    def average_hist_hb(self):
        """
        average histogram number of hb per acceptor, donor

        :return: a numpy array
        """

        hist_avg = np.zeros([self.num_hb_types, 3, 6])
        for ind_time in range(self.num_frames):
            mat = self.compute_hist_hb(ind_time)
            for ind_types in range(self.num_hb_types):
                for ind_g in range(3):
                    hist_avg[ind_types, ind_g, :] += mat[ind_types, ind_g, :]

        return hist_avg/self.num_frames

    def average_hist_hb_peo(self, N_peo=1, peo_size=12):
        """
        average histogram number of hb per acceptor, donor

        :return: a numpy array
        """

        hist_avg = np.zeros([N_peo, 6])
        for ind_time in range(self.num_frames):
            mat = self.compute_hist_hb_per_peo(ind_time, N_peo=N_peo, peo_size=peo_size)
            for ind_types in range(N_peo):
                    hist_avg[ind_types, :] += mat[ind_types, :]

        hist_avg= hist_avg/self.num_frames

        hist_peo_avg = np.zeros([6])
        for ind_types in range(N_peo):
            hist_peo_avg[:] += hist_avg[ind_types, :]

        hist_peo_avg = hist_peo_avg/N_peo

        return hist_peo_avg

    def compute_hist_hb_per_peo(self, index_time, N_peo=1, peo_size=12):
        """
        histogram of hb per acceptor, donor for each frame

        :param index_time: given frame
        :return: a numpy array
        """

        # get tags
        hb_count = self.get_hb_count_per_tag(index_time)
        # get hb histograms for each num_type, at a given time=index_time, for acceptor(0), donor(1) or all(1)
        hb_count = hb_count[self.acceptor_list.index('OP')][0]


        #print(hb_count)
        peo_hb_count = []
        peo_tags = copy.deepcopy(list(self.num_adp[self.acceptor_list.index('OP')][0][0]))
        peo_tags = np.array_split(peo_tags, N_peo)
        #print(peo_tags)
        for peo in peo_tags:
            ind = [list(hb_count[0]).index(oxy) for oxy in peo if oxy in hb_count[0]]
            peo_hb_0 = np.array(hb_count[0])[ind]
            peo_hb_1 = np.array(hb_count[1])[ind]
            peo_hb_count += [(peo_hb_0, peo_hb_1)]

        hist_vals = np.zeros([N_peo, 6])
        for ind_types in range(N_peo):
            for ind_num in range(1, 6):
                hist_vals[ind_types][ind_num] = np.sum(peo_hb_count[ind_types][1] == ind_num)
            # in case there is one molecule with 6 hbs
            num_6 = np.sum(peo_hb_count[ind_types][1] > 5)
            if num_6 > 0:
                print(num_6, 'molecules with 6 or more hydrogen bonds detected')
                hist_vals[ind_types][5] += num_6
            hist_vals[ind_types][0] = peo_size - peo_hb_count[ind_types][1].shape[0]
            hist_vals[ind_types][:] /= peo_size

        return hist_vals

    def compute_hist_hb(self, index_time):
        """
        histogram of hb per acceptor, donor for each frame

        :param index_time: given frame
        :return: a numpy array
        """

        # get tags
        hb_count = self.get_hb_count_per_tag(index_time)
        # get hb histograms for each num_type, at a given time=index_time, for acceptor(0), donor(1) or all(1)
        hist_vals = np.zeros([self.num_hb_types, 3, 6])
        for ind_types in range(self.num_hb_types):
            for ind_g in range(3):
                for ind_num in range(1, 6):
                    hist_vals[ind_types][ind_g][ind_num] = np.sum(hb_count[ind_types][ind_g][1] == ind_num)
                # in case there is one molecule with 6 hbs
                num_6 = np.sum(hb_count[ind_types][ind_g][1] > 5)
                if num_6 > 0:
                    print(num_6, 'molecules with 6 or more hydrogen bonds detected')
                    hist_vals[ind_types][ind_g][5] += num_6
            a_list = self.num_adp[ind_types][0][0]
            d_list = self.num_adp[ind_types][1][0]
            all_list = np.unique(np.concatenate((a_list, d_list)))
            hist_vals[ind_types][0][0] = a_list.shape[0] - hb_count[ind_types][0][1].shape[0]
            hist_vals[ind_types][1][0] = d_list.shape[0] - hb_count[ind_types][1][1].shape[0]
            hist_vals[ind_types][2][0] = all_list.shape[0] - hb_count[ind_types][2][1].shape[0]
            hist_vals[ind_types][0][:] /= a_list.shape[0]
            hist_vals[ind_types][1][:] /= d_list.shape[0]
            hist_vals[ind_types][2][:] /= all_list.shape[0]

        return hist_vals

    def get_hb_count_per_tag(self, index_time):
        """
        tags and number of hb per acceptor, donor for each frame

        :param index_time: given frame
        :return: a numpy array
        """

        lst_vals = [[0] *3  for i in range(self.num_hb_types)]

        for ind_types in range(self.num_hb_types):
            hb_tags = self.hbonds[ind_types][index_time]
            for ind_da in range(2):
                lst_vals[ind_types][ind_da] = np.unique(hb_tags[ind_da], return_counts=True)
            lst_vals[ind_types][2] = np.unique(np.concatenate((hb_tags[0], hb_tags[1])), return_counts=True)

        return lst_vals

    def compute_number_hb(self):
        """
        total number of hydrogen bonds for a given donor, proton and acceptor

        :return: a 2D list containing the number of hydrogen bonds present in each frame for each hb type
        """
        
        return [[self.hbonds[ind_p][ind][0].shape[0] for ind in range(self.num_frames)]
                for ind_p in range(self.num_hb_types)]

    def _get_hbonds(self, donor, proton, acceptor):
        """
        calculates the number of hydrogen bonds given a donor, proton and acceptor
        according to the criteria of Berndt et. al. (1993)

        :param donor: name of the donor molecule
        :param proton: name of the proton
        :param acceptor: name of the acceptor molecule
        :return: a 2D numpy array containing the number of hydrogen bonds present in each frame
        """

        print('Computing Hydrogen Bonds: donor', donor, 'proton ', proton, 'acceptor ', acceptor)

        hb_tags = self.num_frames*[0]
        # distance between donor and proton
        dist_d_p = self.ff.get_bond_r0(donor+'-'+proton)*self.units.length_construction_to_simulation*(1+self.eps)

        # the quantities below cannot change during the trajectory, so we define them outside the loop
        snap = self.sys[self.ini_frame]
        # donor, acceptor and proton names
        a_id = snap.particles.types.index(acceptor)
        d_id = snap.particles.types.index(donor)
        p_id = snap.particles.types.index(proton)
        # identify all possible protons, acceptors and donors
        alla_tags = np.nonzero(snap.particles.typeid == a_id)
        alld_tags = np.nonzero(snap.particles.typeid == d_id)
        allp_tags = np.nonzero(snap.particles.typeid == p_id)

        for ind, ind_frame in enumerate(range(self.ini_frame, self.end_frame)):
            snap = self.sys[ind_frame]
            # get the box and make a box object
            box_sim = snap.configuration.box
            box = BoxFunctions(box_sim)
            # check all acceptor-proton distances that satisfy the hb condition (and remove possible donors)
            pos = snap.particles.position
            p_tags = copy.deepcopy(allp_tags)
            a_tags = copy.deepcopy(alla_tags)
            dist_all = box.compute_all_distances(pos[p_tags], pos[a_tags])
            hb_pa = np.nonzero(np.logical_and(dist_all <= self.r_cut, dist_all > dist_d_p))
            # those are the tags satisfying those conditions
            p_tags = p_tags[0][hb_pa[0]]
            a_tags = a_tags[0][hb_pa[1]]
            # find the tag index for donors candidates
            d_tags = copy.deepcopy(alld_tags)
            hb_pd = np.nonzero(box.compute_all_distances(pos[p_tags], pos[d_tags]) < dist_d_p)
            d_tags = d_tags[0][hb_pd[1]]
            if len(a_tags) != len(d_tags) or len(a_tags) != len(p_tags):
                raise ValueError('inconsistency found in hb calculation')
            # we need to keep only those that satisfy the angle constrain
            dist1 = box.compute_distances(pos[d_tags], pos[a_tags])
            dist2 = box.compute_distances(pos[d_tags], pos[p_tags])
            dist3 = box.compute_distances(pos[a_tags], pos[p_tags])
            cos_val = -(dist1*dist1-dist2*dist2-dist3*dist3)/(2*dist2*dist3)
            ind_hb = np.logical_and(cos_val >= self.min_cos, cos_val <= self.max_cos)
            p_tags = p_tags[ind_hb]
            a_tags = a_tags[ind_hb]
            d_tags = d_tags[ind_hb]
            # those are the tags (acceptor, donor, proton) defining the hydrogen bonds
            hb_tags[ind] = (a_tags, d_tags, p_tags)

        return hb_tags, (alla_tags, alld_tags, allp_tags)
