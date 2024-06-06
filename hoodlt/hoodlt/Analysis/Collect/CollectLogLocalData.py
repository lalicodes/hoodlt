"""
:module: CollectLogLocalData
:Platform: Windows, Unix
:synopsis: Gets Per particle data from log files to be used in simulation analysis and calculations

.. moduleauthor:: Alex Travesset, August 2022
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, January 2024
..                        rewrote the code to be consistent with the changes in analyze
"""
import gsd.hoomd
import numpy as np
import glob
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader


class CollectLogLocalData:
    """
    This class gets data from the log files given in the constructor.
    """
    def __init__(self, file_name, ff, num_frames=None, end_val=None):
        """
        Creates a CollectLogData object. 

        :param file_name: name of .gsd file to get the data from, without the _local_number_.gsd extension
        :param ff: name of the force field used to build and run the simulation
        :param num_frames: number of frames to consider
        :param end_val: end value
        """

        self.file_name = file_name
        f_name = self.file_name + '_local_' + '.gsd'

        self.file_traj = gsd.hoomd.open(f_name, 'r')

        self.units = ForceFieldReader(ff).get_units()

        self.num_frames = len(self.file_traj)

        if end_val is None:
            self.end_val = self.num_frames
        else:
            self.end_val = end_val

        if num_frames is not None:
            self.num_frames = num_frames

        self.ini_val = self.end_val - self.num_frames

        if self.ini_val < 0:
            raise ValueError('the number of frames is larger than the available number')

        self.quant_hoomd = list(self.file_traj[self.ini_val].log.keys())
        self.num_quantities = len(self.quant_hoomd)
        self.shape = list(self.file_traj[self.ini_val].log[self.quant_hoomd[-1]].shape)
        # quantities coming up from hoomd are too verbose, so we use simpler names and store them in quant

        dict_names = {'md/': '', 'particles/': '', '/': '_'}

        self.quant = len(self.quant_hoomd)*['c']
        for ind, val in enumerate(self.quant_hoomd):
            self.quant[ind] = val
            for st_original, st_replace in dict_names.items():
                self.quant[ind] = self.quant[ind].replace(st_original, st_replace)

    def loggables(self):
        """
        Returns all the available loggables
        """

        return self.quant

    def get_quantity(self, quantity, nc_entity='all'):
        """
        Gets every data point for the given quantity in every one of the list of files

        :param quantity: the quantity to get from the list of loggabbles
        :param nc_entity: a nc entity
        :return: a list of numpy array with dimensions [number of files, number of data points in a file]
        """

        indx = self.quant.index(quantity)

        if nc_entity == 'all':
            n_i = 0
            n_f = self.file_traj[self.ini_val].log[self.quant_hoomd[indx]].shape[0]
        else:
            n_i = nc_entity['tag-center']
            n_f = nc_entity['tag-end']+1

        list_res = []
        for ind in range(self.ini_val, self.end_val):
            val = self.file_traj[ind].log[self.quant_hoomd[indx]][n_i:n_f]
            list_res.append(val)

        return np.array(list_res)
