"""
:module: CollectGeneralData
:Platform: Windows, Unix
:synopsis: Class to collect data for analysis from sources other than forcefields, gsds, log files, or hist files

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
"""

import numpy as np


class CollectGeneralData:
    """
    This class collects data from any files that aren't gsd, log, forcefield, or hist files. It will return values
    as they appear in the files, and the user can get those values scaled by a constant, as usual
    """
    def __init__(self, filename):
        """

        :param filename: the name of the file to get the data from, with all file extensions
        """

        self.filename = filename

    def get_columns_from_text_file(self):
        """
        Gets all the data from all the columns in a text file

        :return: a 2D numpy array with dimensions [number of columns, number of data points per column]
        """

        data = np.genfromtxt(self.filename)

        result = []
        for i in range(len(data[0])):
            result.append(data[:,i])

        return np.array(result)
