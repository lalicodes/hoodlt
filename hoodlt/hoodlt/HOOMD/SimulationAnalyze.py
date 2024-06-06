"""
:module: SimulationAnalyze
:platform: Unix, Windows
:synopsis: Class which computes local and other more complex quantities needed to analyze results

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, August2022
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, January 2024
..                        adapted code to compute forces
"""

import hoomd
from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation


class SimulationAnalyze(HoomdSimulation):
    """
    Class designed to compute local and other quantities
    """

    def __init__(self, sysfile, rcut, ff_name, cluster=False, temp_in_kelvin=None, seed=5, dict_params={},
                 dict_nlist=None, list_electrostatic=None, special_rcut={}, translation=True,
                 rotation=True, walls_dictionary=None, drift_removal_dict=None):
        """
        Simulation recording bond fluctuations necessary for wham and free energy

        :param sysfile: systemfile
        :param rcut: the rcut for the system (i.e. each nonbonded pair will have cutoff rcut*sigma)
        :param ff_name: name of the force field,
        :param cluster: cluster option
        :param temp_in_kelvin: temperature in kelvin units
        :param seed: seed
        :param dict_nlist: option to override the neighborlist. Should be either 'tree' 'cell' or 'Stencil'
        :param dict_params: integrator dictionary,
        :param special_rcut: dictionary for specifying nonbonded rcuts which are not the same as the rcut argument.
        :param translation: whether to integrate rigid bodies
        :param rotation: whether to integrate rotational degree
        :param walls_dictionary: Given as a dictionary to make walls
        :param drift_removal_dict: sets drift removal to zero
        """

        # start the new class
        super(SimulationAnalyze, self).__init__(sysfile, rcut, ff_name, cluster=cluster, temp_in_kelvin=temp_in_kelvin,
                                                seed=seed, dict_params=dict_params, write_data=True,
                                                dict_nlist=dict_nlist, list_electrostatic=list_electrostatic,
                                                special_rcut=special_rcut, translation=translation, rotation=rotation,
                                                walls_dictionary=walls_dictionary,
                                                drift_removal_dict=drift_removal_dict)

    def run_hoomd_analyze(self, sim_param, dt=None, log_partial_energies=False, quants=None):
        """
        runs the simulation

        :param sim_param: SimParameters object
        :param dt: time step
        :param log_partial_energies: whether to log the partial contributions to the potential energy.
        :param quants: this is a dictionary, key is the name of the object, value are the quantities to log, quantities
            are precisely defined by hoomd, for example: energy, energies, forces, torques, etc.. the key is the name of
            the object. For example, pair_lj, bondharmonic, etc.. as defined in hoomd simulation from hoodlt.
        """

        if quants is None:
            log_quantities = ['energies']
        else:
            log_quantities = quants

        file_name = self.sysfile + '_local_'+'.gsd'
        l_logs = hoomd.logging.Logger()
        trig = hoomd.trigger.Periodic(sim_param.steps_log)

        for k_pots, pots in self.potentials.items():
            if k_pots in log_quantities:
                l_logs.add(pots, quantities=log_quantities[k_pots])

        writer_local = hoomd.write.GSD(filename=file_name, trigger=trig, mode='ab', filter=hoomd.filter.All())
        writer_local.logger = l_logs

        self.gsd_writers.append(writer_local)

        self.run_hoomd(sim_param, dt=dt, log_partial_energies=log_partial_energies)
