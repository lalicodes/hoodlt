"""
:module: CollectLogData
:Platform: Windows, Unix
:synopsis: Gets data from log files to be used in simulation analysis and calculations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - eliminated explicit references to units classes
..                Alex Travesset <trvsst@ameslab.gov>, April 2022
..                  - made the class compatible with hoomd v3 and above
..                Alex Travesset <trvsst@ameslab.gov>, July 2022
..                  - Added more rigorous error controls to the files
..
"""
import gsd.hoomd
import numpy as np
from scipy.stats import sem
from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader


class CollectLogData:
    """
    This class gets data from the log files given in the constructor.
    """
    def __init__(self, file_names, ff, num_frames=None, typ='therm', end_val=None, quant_only=(-1, -1)):
        """
        Creates a CollectLogData object. 

        :param file_names: list of the names of .gsd files to get the data from, without the _therm.gsd extension
        :param ff: name of the force field used to build and run the simulation
        :param num_frames: number of frames to consider
        :param: file type, either therm (thermodynamic quantities) or sim (simulation status quantities)
        :param end_val: end value
        :param quant_only: return just two quantities (in the order stored) is useful to just return a subset of data
        """

        self.file_names = [fn+'_'+typ for fn in file_names]
        self.file_traj = [gsd.hoomd.open(fn+'.gsd', 'r') for fn in self.file_names]
        self.units = ForceFieldReader(ff).get_units()

        self.num_files = len(file_names)
        self.num_data = [len(self.file_traj[ind]) for ind in range(self.num_files)]

        min_data = min(self.num_data)
        if end_val is None:
            self.end_val = self.num_data
        else:
            self.end_val = [end_val]*self.num_files

        if num_frames is None:
            num_frames = min_data
        if num_frames > min_data:
            raise ValueError('Incorrect number of frames. The maximum number of frames possible is %d' % min_data)

        self.num_frames = num_frames
        self.ini_val = [e_val - self.num_frames for e_val in self.end_val]

        if quant_only == (-1, -1):
            speed_up = False
        else:
            speed_up = True

        self.quant_hoomd = list(self.file_traj[0][self.ini_val[0]].log.keys())
        self.num_quantities = len(self.quant_hoomd)
        # quantities coming up from hoomd are too verbose, so we use simpler names and store them in quant
        if typ == 'therm':
            dict_names = {'md/': '', 'compute/ThermodynamicQuantities/': '', '/': '_'}
        else:
            dict_names = {'md/': '', '/': '_'}

        self.quant = len(self.quant_hoomd)*['c']
        for ind, val in enumerate(self.quant_hoomd):
            self.quant[ind] = val
            for st_original, st_replace in dict_names.items():
                self.quant[ind] = self.quant[ind].replace(st_original, st_replace)

        # collect the data as specified and place it in a numpy array
        log_data = np.zeros([self.num_files, self.num_frames, self.num_quantities])

        for ind_f, file_name in enumerate(self.file_names):
            i_v = self.ini_val[ind_f]
            i_f = self.end_val[ind_f]
            for ind_q, quant in enumerate(self.quant_hoomd):
                if speed_up and not ind_q in quant_only:
                    pass
                else:
                    for ind_d, ind_val in enumerate(range(i_v, i_f)):
                        log_data[ind_f, ind_d, ind_q] = self.file_traj[ind_f][ind_val].log[quant]

        self.log_data = log_data

    def get_quantity(self, quantity):
        """
        Gets every data point for the given quantity in every one of the list of files

        :param quantity: the quantity to get
        :return: a 2D numpy array with dimensions [number of files, number of data points in a file]
        """

        ind_quantity = self.quant.index(quantity)

        return self.log_data[:, :, ind_quantity]

    def get_average(self, quantity):
        """
        Gets the average of the given quantity in every one of the list of files

        :param quantity: the quantity to get
        :return: a tuple of 1D numpy array with dimensions [number of files] with average, sem and std
        """

        ind_quantity = self.quant.index(quantity)
        average = np.average(self.log_data[:, :, ind_quantity], axis=1)
        st_dev = np.std(self.log_data[:, :, ind_quantity], axis=1)
        st_err = sem(self.log_data[:, :, ind_quantity], axis=1)

        return [(average[ind], st_dev[ind], st_err[ind]) for ind in range(average.shape[0])]
