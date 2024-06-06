"""
:module: SimParameters
:platform: Unix, Windows
:synopsis: Implements mpi mean force calculation

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
"""


class SimParameters(object):
    """
    Simple class containing the parameters of a simulation
    """

    def __init__(self, steps_sim, steps_log, quantities_log, steps_write, steps_hist=None, num_procs=1):
        """The constructor

        :param steps_sim: number of steps in the simulation
        :param steps_log: number of steps after which observables are measured
        :param quantities_log: list of quantities to be logged
        :param steps_write: number of steps after which configurations are saved
        :param steps_hist: number of steps after which to output histogramming
        :param num_procs: number of processors
        """

        self.steps_sim = steps_sim
        self.steps_log = steps_log
        self.quantities_log = quantities_log
        self.steps_write = steps_write
        self.steps_hist = steps_hist
        self.num_procs = num_procs
