"""
:module: CollectHistData
:Platform: Windows, Unix
:synopsis: Gets data from a .hist file to be used in simulation analysis and calculations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - eliminated explicit references to units classes
..                  - Updated Documentation
..                Alex Travesset <trvsst@ameslab.gov>, April 2022
..                  - Made it consistent with hoomd v3
"""

import numpy as np
import json
import gsd.hoomd
import hoodlt.Data.Modelconfigurations.Saver as Sv


class CollectHistData:
    """
    This class collects data from hist files, and returns data in simulation units (not dimensionless)
    """
    def __init__(self, file_names, num_frames=None, end_val=None):
        """

        :param file_names: list of hist files to get data from, without the '_hist.gsd' extension
        :param num_frames: number of frames to consider
        :param end_val: Instead of taking the last num_frames in the series, skip the last end_val frames.
        """

        self.file_names = file_names
        file_traj = [gsd.hoomd.open(fn + '_hist.gsd', 'r') for fn in self.file_names]
        # number of files
        self.num_files = len(self.file_names)
        # number of frames on each file
        self.num_data = [len(file_traj[ind]) for ind in range(len(self.file_names))]

        # establish  the initial and final frames to use
        min_data = min(self.num_data)
        self.end_val = self.num_data
        if isinstance(end_val, int):
            self.end_val = [end_val]*self.num_files
        elif type(end_val) is list:
            self.end_val = end_val

        if num_frames is None:
            num_frames = min_data
        self.num_frames = num_frames
        self.ini_val = [e_val - self.num_frames for e_val in self.end_val]

        # check that the number of frames requested is consistent with the frames available
        if min_data < num_frames:
            if any(ind < 0 for ind in self.ini_val):
                print('number of frames ', min_data, 'number of frames to average ', num_frames)
                print('initial values')
                print(self.ini_val)
                raise ValueError('The number of frames exceeds the frames available')

        # make a dictionary of the first file and extract the parameter values
        params = ['types', 'temperature', 'bond histogram']
        self.dict_bond_data = {}
        with open(self.file_names[0] + '_bonds.json') as fp:
            dict_bond_data = json.load(fp)
            for key in params:
                self.dict_bond_data[key] = dict_bond_data[key]

        # define the spring distance and lamda values, there are as many distances as bond types
        self.bond_r0 = np.zeros([self.num_files, len(self.dict_bond_data['types'])])
        self.lamda = np.zeros([self.num_files, len(self.dict_bond_data['types'])])

        # read spring distances and check that the parameters are the same for all files
        for ind_f, fl in enumerate(self.file_names):
            with open(fl + '_bonds.json') as fp:
                dict_bond_data = json.load(fp)
                self.bond_r0[ind_f, :] = dict_bond_data['dist'][:]
                self.lamda[ind_f, :] = dict_bond_data['lamda'][:]
                for val in params:
                    if dict_bond_data[val] == self.dict_bond_data[val]:
                        pass
                    else:
                        print(val, dict_bond_data[val], self.dict_bond_data[val])
                        raise ValueError('json files contain inconsistent parameters')
    
        # define numpy array and get the bond fluctuations
        self.num_types = len(dict_bond_data['types'])
        if dict_bond_data['bond histogram'] == 'average':
            func_type = 'average'
            self.dim_bonds = None

        elif dict_bond_data['bond histogram'] == 'all':
            func_type = 'all'
            # read the number of bonds
            conf = Sv.load_config(self.file_names[0])
            self.dim_bonds = np.zeros(self.num_types, dtype=int)
            for ind in range(self.num_types):
                self.dim_bonds[ind] = len(conf.bonds[ind])
            self.num_types = np.sum(self.dim_bonds)

        self.ts = np.zeros([self.num_files, self.num_frames])
        self.bond_data = np.zeros([self.num_files, self.num_frames, self.num_types])
        for ind_f, fl in enumerate(file_traj):
            ini = self.ini_val[ind_f]
            end = self.end_val[ind_f]
            for ind_t, frame in enumerate(fl[ini:end]):
                self.ts[ind_f, ind_t] = frame.log['Simulation/timestep']
                for ind_l in range(self.num_types):
                    self.bond_data[ind_f, ind_t, ind_l] = frame.log['ParticleDistance/'+func_type][ind_l]

    def get_fluctuations(self, num):
        """
        obtains changes in bond lengths

        :param num: index of the file for which the bonds are given
        :return : numpy array for the given bonds [bonds at a given time, number of bonds]
        """

        return self.bond_r0[num], self.bond_data[num]

    def get_timestep(self):
        """
        provides the timesteps at which the bonds are taken

        :return: numpy array with the time steps
        """

        return self.ts
