"""
:module: SimulationWithBonds
:platform: Unix, Windows
:synopsis: Class which can initialize the CTR-CTR bonds in a simulation

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - can now use any type of hoomd 
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - units are now taken from ff reader instead of as an argument
..                  - Added documentation
..                  - Nonbonded rcuts are now set pairwise, instead of using a system wide value
..                  - Certain nonbonded interactions can now have their energies scaled
..                Alex Travesset <trvsst@ameslab.gov>, April 2022
..                  - made it consistent with hoomd v3
"""
import shutil
import json
import hoomd
from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
from hoodlt.HOOMD.SimulationActions import ParticleDistance


class SimulationWithBonds(HoomdSimulation):
    """
    Class designed run simulations while also being able to manipulate the CTR-CTR bonds in the system
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
        :param nlist: option to override the neighborlist. Should be either 'tree' 'cell' or 'Stencil'
        :param dict_params: simulation parameters for the NVT integrator as a dictionary
        :param special_rcut: dictionary for specifying nonbonded rcuts which are not the same as the rcut argument.
        :param translation: whether to integrate rigid bodies
        :param rotation: whether to integrate rotational degree
        :param walls_dictionary: Given as a dictionary to make walls
        :param drift_removal: sets drift removal to zero
        """

        # start the new class
        super(SimulationWithBonds, self).__init__(sysfile, rcut, ff_name, cluster=cluster,
                                                  temp_in_kelvin=temp_in_kelvin, seed=seed, dict_params=dict_params,
                                                  write_data=True, dict_nlist=dict_nlist,
                                                  list_electrostatic=list_electrostatic, special_rcut=special_rcut,
                                                  translation=translation, rotation=rotation,
                                                  walls_dictionary=walls_dictionary, drift_removal_dict=drift_removal_dict)

        self.local_snapshot = 'cpu_local_snapshot'
        if cluster:
            self.local_snapshot = 'gpu_local_snapshot'

        # bond data json file must be generated from ConfigurationBuilder
        with open(sysfile + '_bonds.json') as fp:
            self.dict_bond_data = json.load(fp)

        # checking that there is no mismatch between json and restart (gsd) files
        if 'restart_name' not in self.dict_bond_data:
            if self.dict_bond_data['name'] != sysfile:
                raise ValueError('file gsd name does not match the one contained in .json')

        if 'restart_name' in self.dict_bond_data:
            if self.dict_bond_data['restart_name'] != sysfile:
                raise ValueError('file gsd name does not match the one contained in .json')
            self._init_bonds()

    def initialize(self, add_to_file_name, bond_name=None, rinit=None, lamda=None, log_hist=None):
        """
        initializes the parameters for the bond type bond_name

        :param add_to_file_name: string to add to the file name
        :param bond_name: list of bond names
        :param rinit: list of r0 bond distances
        :param lamda: list of dimensionless scaling factor for bond strength, defaults to 1
        :param log_hist: None (do not log hist), all (log each individual bond), average
        """

        if bond_name is None:
            bond_name = self.dict_bond_data['types']

        if rinit is None:
            rinit = self.dict_bond_data['dist']

        for bnd_name in bond_name:
            if bnd_name not in self.dict_bond_data['types']:
                raise ValueError('the list of bonds provided does not match')

        if 'restart_name' in self.dict_bond_data:
            print('bonds can only be initialized once. Skipping reinitialization ')
        else:
            self.dict_bond_data['dist'] = rinit
            if lamda is None:
                lamda = len(rinit)*[1.0]
            self.dict_bond_data['lamda'] = lamda

            # make sure that from now on, the file cannot be modified and give the new name
            self.dict_bond_data['restart_name'] = self.sysfile + '_' + add_to_file_name

            # type of bonds to log
            accepted_hist_names = ['none', 'all', 'average']

            if log_hist not in accepted_hist_names:
                log_hist = accepted_hist_names[0]

            self.dict_bond_data['bond histogram'] = log_hist
            self.dict_bond_data['temperature'] = self.realtemp

            new_name = self.dict_bond_data['restart_name'] + '_bonds.json'
            with open(new_name, 'w') as fp:
                json.dump(self.dict_bond_data, fp)

            self.file_write = self.dict_bond_data['restart_name'] + '_write.gsd'
            self.file_restart = self.dict_bond_data['restart_name'] + '_restart.gsd'
            self.file_sim = self.dict_bond_data['restart_name'] + '_sim.gsd'
            self.file_thermo = self.dict_bond_data['restart_name'] + '_therm.gsd'
            f_pickle = self.file_pickle
            self.file_pickle = self.dict_bond_data['restart_name'] + '.pickle'
            shutil.copyfile(f_pickle, self.file_pickle)
            self._init_bonds()

    def _init_bonds(self):
        """
        initializes bonds internally
        """

        for ind, bond_name in enumerate(self.dict_bond_data['types']):
            r1 = self.dict_bond_data['dist'][ind]
            k_ff = self.ff_reader.get_potentials_params('bond', 'harmonic', bond_name)
            k_temp = k_ff['k'] * self.dict_bond_data['lamda'][ind]
            self.potentials['bondharmonic'].params[bond_name] = dict(r0=r1, k=k_temp)

    def run_hoomd_hist(self, sim_param, dt=None, log_partial_energies=False, additional_filters=None):
        """
        runs the simulation

        :param sim_param: SimParameters object
        :param dt: time step
        :param log_partial_energies: whether to log the partial contributions to the potential energy.
        :param additional_filters: list of additional filters
        """

        if 'restart_name' not in self.dict_bond_data:
            raise ValueError('System needs to be initialized, use member function initialize()')

        real_temperature = self.dict_bond_data['temperature']
        if real_temperature != self.realtemp:
            raise ValueError('Use standard protocols to change temperature')

        self.file_dict['hist_params'] = self.dict_bond_data['restart_name'] + '_bonds.json'

        sim_hist_steps = sim_param.steps_hist
        if sim_hist_steps is None:
            sim_hist_steps = sim_param.steps_log

        log_hist = self.dict_bond_data['bond histogram']
        if log_hist != 'none':
            self.file_dict['hist_data'] = self.dict_bond_data['restart_name'] + '_hist.gsd'
            logger4 = hoomd.logging.Logger()
            logger4.add(self.sim, quantities=['timestep'])
            trig = hoomd.trigger.Periodic(sim_hist_steps)
            nm = self.file_dict['hist_data']
            gsd_writer4 = hoomd.write.GSD(filename=nm, trigger=trig, mode='ab', filter=hoomd.filter.Null())
            gsd_writer4.logger = logger4

            bond_params = self.dict_bond_data['types']
            p_dist = ParticleDistance(self.sim, bond_params, self.local_snapshot)
            logger4[('ParticleDistance', log_hist)] = (p_dist, log_hist, 'sequence')
        
            self.gsd_writers.append(gsd_writer4)

        self.run_hoomd(sim_param, dt=dt, log_partial_energies=log_partial_energies, additional_filters=additional_filters)
