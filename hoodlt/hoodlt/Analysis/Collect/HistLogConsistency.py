"""
:module: HistLogConsistency
:Platform: Windows, Unix
:synopsis: Checks consistency between log and hist files

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, July 2022
.. history:
"""

import numpy as np


def consistency_therm_and_hist_files(coll_log, coll_hist):
    """
    check that the energies and histograms are taken at the same time steps
    """

    ts_therm = coll_log.get_quantity('Simulation_timestep')
    ts_hist = coll_hist.get_timestep()

    if not np.array_equiv(ts_therm, ts_hist):
        raise ValueError('history and energy at taken at different time steps')