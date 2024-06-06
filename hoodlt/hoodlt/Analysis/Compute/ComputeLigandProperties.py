"""
:module: ComputeLigandProperties
:Platform: Windows, Unix
:synopsis: Computes Properties of the Ligands in a simulation

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - class is now constructed using an explicit trajectory object
..                  - eliminated explicit references to units classes
..                  - updated documentation
..                  - Moved old code (written by Curt Waltmann) to fit new analysis scheme
..
..                Alex Travesset <trvsst@ameslab.gov>, August 2022
..                  - Revamped the class to include further analysis
..
"""

import numpy as np


class ComputeLigandProperties():
    """
    This class compute properties of ligands within nanocrystals
    """

    def __init__(self, nc):
        """

        :param nc: Functionalized particle object
        """

        self.nc = nc

    def compare_ligands(self, nc2):
        """
        returns the average positions
        """

        mat_diff = []

        for lig1, lig2 in zip(self.nc.ligands, nc2.ligands):
            c0 = lig1.position[0]
            c1 = lig2.position[0]
            diff = np.subtract(lig1.position - c0, lig2.position - c1)
            mat_diff.append(diff)

        return mat_diff
