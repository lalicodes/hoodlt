"""
:module: HoomdSimulation
:platform: Unix
:synopsis: Defines the class to prepare a hoomd simulation

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2017
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, July 2017
..                        simplified code
..                Tommy Waltmann <tomwalt@iastate.edu>, April 2019
..                  - Rewrote constructor to interface with the new forcefield reader classes
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - units are now taken from ff reader instead of as an argument
..                  - added documentation
..                  - added impropers
..                  - Nonbonded rcuts are now set pairwise, instead of using a system wide value
..                  - Certain nonbonded interactions can now have their energies scaled
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - can now use any type of hoomd integrator
..                Jacob Austin <jaustin2@iastate.edu> April 2020
..                  - added special coulomb and lj functionality
..                Xun Zha <xzha@iastate.edu> Dec 2020
..                  - changed default value of the parameter dscale from False to 0.0116434
..                Alex Travesset <trvsst@ameslab.gov> June 2021
..                  - changed the timestep so that it is expressed in simulation units
..                  - ensured that self.units.dt is the actual timestep being used
..                Jianshe Xia <xiajs6075@iccas.ac.cn> July 2021
..                  - added the table force field for bond, angle and dihedral
..                Jonas Hallstrom <jlhallst@asu.edu> July 2021
..                  - added rotation parameter to control anisotropic integration
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                  - upgraded the class to HOOMD v3 and beyond
..                  - significantly simplified code and eliminated unnecessary member functions
..                Elizabeth Macias <emacias@iastate.edu> August 2022
..                  - added implemenation for reading, setting, and saving degrees-of-freedom parameters
..                Elizabeth Macias <emacias@iastate.edu> September 2022
..                  - implemented center of mass drift removal
..                Alex Travesset <trvsst@ameslab.gov> September 2023
..                  - made it compatible with HOOMD v4 and beyond
..                Alex Travesset <trvsst@ameslab.gov> November 2023
..                  - Simplified the handling of thermostats, thermostat file is now a json file
"""
import copy
import os
import os.path
import json
import copy
import hoomd
import hoodlt.HOOMD.CoefficientSetter as cs
from hoodlt.HOOMD.GsdParser import GsdParser
from hoodlt.Data.Forcefield.ScaledForceFieldReader import ScaledForceFieldReader


class HoomdSimulation(object):
    """
    Class which creates simulation contexts and runs simulations using HOOMD
    """

    def __init__(self, sysfile, rcut, ff_name, cluster=False, temp_in_kelvin=None, seed=5, dict_params={},
                 write_data=True,  dict_nlist=None, list_electrostatic=None,
                 special_rcut={}, translation=True, rotation=True, walls_dictionary=None, drift_removal_dict=None):
        """
        Creates the context for running simulations

        :param sysfile: systemfile
        :param rcut: the rcut for the system (i.e. each nonbonded pair will have cutoff rcut*sigma)
        :param ff_name: name of the force field
        :param cluster: cluster option activates the GPU
        :param temp_in_kelvin: temperature in  kelvin
        :param dict_params: simulation_parameters for the NVT integrator as a dictionary
        :param write_data: save data in files
        :param dict_nlist: dictionary option to override the neighborlist
        :param list_electrostatic: parameter list of for pppm electrostatic: None or [resolution, order, rcut, alpha]
        :param special_rcut: cutoff for special interactions, if specified needs to be a dictionary
        :param translation: whether to integrate rigid bodies to  altogether
        :param rotation: whether to integrate rotational degrees of freedom
        :param walls_dictionary: Given as a dictionary to make walls
        :param drift_removal_dict: sets drift removal to zero
        """

        # definitions
        self.sysfile = sysfile
        self.rot = rotation
        self.seed = seed
        self.write_data = write_data

        # define the default dictionary defining the neighborlist
        if dict_nlist is None:
            dict_nlist = {'nlist': 'Cell', 'params': {'buffer': 0.4}}

        # keep track of the potentials that are defined
        self.potentials = {}

        # globalizing drift removal dictionary
        self.drift_removal_dict = drift_removal_dict

        # define the device
        if isinstance(cluster, dict):
            com = cluster['com']
            print('partition ', com.partition, ' with rank ', com.rank)
            self.sim_dev = hoomd.device.GPU(communicator=com)
        elif cluster:
            self.sim_dev = hoomd.device.GPU()
        else:
            self.sim_dev = hoomd.device.CPU()

        # make a force field reader
        self.ff_reader = ScaledForceFieldReader(ff_name)
        self.units = self.ff_reader.get_units()

        # read the files and run only if there is a restart file
        self.sim = hoomd.Simulation(self.sim_dev, seed=self.seed)
        self.sim.create_state_from_gsd(filename=self.sysfile + '_restart.gsd')
        self.snap = self.sim.state.get_snapshot()
        self.gsd_parser = GsdParser(self.snap)

        # collect degrees of freedom data from file
        self.dict_thermostat_dof = {}
        if os.path.isfile(self.sysfile + '_restart_previous.gsd'):
            with open(self.sysfile + '_restart_thermostat.json') as tp:
                self.dict_thermostat_dof = json.load(tp)

        # define rigid bodies, if necessary
        if self.gsd_parser.has_rigid_bodies():
            self.rigid_body = hoomd.md.constrain.Rigid()
            cs.set_rigid_bodies(self.rigid_body, self.snap, self.sysfile)

        # next need to group particles for integration
        self.particles = hoomd.filter.All()
        self.centers = hoomd.filter.Rigid(('center',))
        self.nonrigid = hoomd.filter.Rigid(('free',))
        self.part_dof = hoomd.filter.Union(self.centers, self.nonrigid)

        # set the neighborlist, whether tree or cell
        nlist = dict_nlist['nlist']
        # exclusions are determined by the force field
        dict_nlist['params']['exclusions'] = self.ff_reader.get_exclusions()
        if nlist == 'Tree' or nlist == 'Cell' or nlist == 'Stencil':
            self.nl = getattr(hoomd.md.nlist, nlist)(**dict_nlist['params'])
        else:
            raise ValueError('nlist should be either Tree, Cell or Stencil')

        # initialize all attributes
        attribs = self.ff_reader.get_list_attributes()
        # particles have been initialized in snap, remove groups
        attribs.remove('groups')
        # set up electrostatic
        if self.gsd_parser.has_charges():
            dict_e = cs.set_electrostatic(list_electrostatic, rcut)
            self.pppm = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(self.nl, **dict_e)
            self.potentials['pppm0'] = self.pppm[0]
            self.potentials['pppm1'] = self.pppm[1]

        # set up non-bonded terms
        for pair in self.ff_reader.get_potentials('nonbonded'):
            pair_name = 'pair_' + pair  # for example, pair_lj
            # get the force field parameters
            p_obj = cs.set_nonbonded(hoomd.md.pair, self.nl, pair, self.snap, self.ff_reader, rcut)
            setattr(self, pair_name, p_obj)  # for example, self.pair_lj = hoomd.md.pair.LJ(nlist)
            self.potentials[pair_name] = p_obj
        # non-bonded pairs have been initialized remove them
        attribs.remove('nonbonded')

        # iterate over all bonded terms
        list_bonded = self.ff_reader.get_list_attributes_without_cutoff()
        for attr in list_bonded:
            state = getattr(self.snap, attr+'s')
            # loop if there is at least a non-zero element
            if getattr(state, 'N') > 0:
                hm_obj = getattr(hoomd.md, attr)  # for example hm_obj = hoomd.md.bond
                for potential in self.ff_reader.get_potentials(attr):
                    # get the force field parameters
                    p_obj = cs.set_bonded(attr, hm_obj, state, potential, self.ff_reader)
                    self.potentials[attr+potential] = p_obj

        # remove all bonded terms, what remains are terms with cutoff (typically special pairs)
        attribs = list(set(attribs)-set(list_bonded))
        for attr in attribs:
            if attr == 'special_pair':
                state = getattr(self.snap, 'pairs')
            else:
                state = getattr(self.snap, attr+'s')
            if getattr(state, 'N') > 0:
                hm_obj = getattr(hoomd.md, attr)  # for example hm_obj1 = hoomd.md.special_pair
                for potential in self.ff_reader.get_potentials(attr):
                    s_cut = rcut
                    if potential in special_rcut:
                        s_cut = special_rcut[potential]
                    p_obj = cs.set_special(attr, hm_obj, state, potential, self.ff_reader, s_cut)
                    self.potentials[attr+potential] = p_obj

        # add walls if specified
        if walls_dictionary is not None:
            # make the walls
            list_walls = walls_dictionary['walls']
            walls_included = []
            for dict_w in list_walls:
                for key, params in dict_w.items():
                    walls_included.append(getattr(hoomd.wall, key)(**params))
            # define the potential
            lj_wall = hoomd.md.external.wall.LJ(walls=walls_included)
            lj_wall.params[self.snap.particles.types] = {'sigma': 0.0, 'epsilon': 0.0, 'r_cut': 0.0}
            list_params = walls_dictionary['params']
            for dict_p in list_params:
                for key, params in dict_p.items():
                    lj_wall.params[key] = params
            # add the potential
            self.potentials['wall'] = lj_wall

        # define the integration group
        groups = [self.particles, self.part_dof, self.nonrigid]
        if translation:
            self.int_group = groups[self.gsd_parser.integration_group()]
        else:
            self.int_group = groups[2]
        # define the integrator
        self.integrator = hoomd.md.Integrator(dt=self.units.dt)

        if self.gsd_parser.has_rigid_bodies():
            self.integrator.integrate_rotational_dof = self.rot
            self.integrator.rigid = self.rigid_body

        if self.snap.constraints.N > 0:
            dist_constrain = hoomd.md.constrain.Distance()
            self.integrator.constraints.append(dist_constrain)

        # add the chosen integrator
        if len(self.dict_thermostat_dof) > 0:
            # if the temperature is set in Kelvin that superseeds any other choice
            if temp_in_kelvin is not None:
                self.realtemp = self.units.simulation_temp(temp_in_kelvin)
            self.change_integration_method(self.dict_thermostat_dof, read_from_file=True)
        # or set the default
        else:
            # set the nvt integrator parameters
            # set temperature
            self.realtemp = self.units.sim_temp_default
            if 'kT' in dict_params:
                self.realtemp = dict_params['kT']
            # choose value for tau
            if 'tau' in dict_params:
                self.tau = dict_params['tau']
            else:
                tau_recommended = 100*self.units.dt
                self.tau = tau_recommended
            # specifying temperature in Kelvin overrules all the other cases
            if temp_in_kelvin is not None:
                self.realtemp = self.units.simulation_temp(temp_in_kelvin)
            self.dict_thermostat_dof['name'] = 'ConstantVolume'
            self.dict_thermostat_dof['thermostat'] = 'MTTK'
            self.dict_thermostat_dof['params'] = {}
            self.dict_thermostat_dof['th_params'] = {'kT': self.realtemp, 'tau': self.tau}
            mttk = hoomd.md.methods.thermostats.MTTK(kT=self.realtemp, tau=self.tau)
            self.method = hoomd.md.methods.ConstantVolume(filter=self.int_group, thermostat=mttk)

        # set up files to store configurations, final state and log simulation and thermodynamics
        self.file_write = self.sysfile + '_write.gsd'
        self.file_restart = self.sysfile + '_restart.gsd'
        self.file_sim = self.sysfile + '_sim.gsd'
        self.file_thermo = self.sysfile + '_therm.gsd'
        self.file_pickle = self.sysfile + '.pickle'
        self.file_original_restart = None
        self.gsd_writers = []
        self.file_dict = {}
        # add integrator
        self.integrator.methods.append(self.method)
        self.sim.operations.integrator = self.integrator

        # add drift removal
        if self.drift_removal_dict:
            trig_list = [getattr(hoomd.trigger, trig)(int(self.drift_removal_dict['triggers'][trig])) for trig in self.drift_removal_dict['triggers']]
            if self.drift_removal_dict['logic'] == None:
                drift_trigger = trig_list[0]
            else:
                drift_trigger = getattr(hoomd.trigger, self.drift_removal_dict['logic'])(trig_list)
            self.sim.operations.updaters.append(hoomd.md.update.ZeroMomentum(drift_trigger))

    def change_integration_method(self, method_dict, read_from_file=False):
        """
        Changes the type of integration method used in the simulation

        :param method_dict: method object
        :param read_from_file: whether this is a thermostat read from file or set up directly from method_dict
        :return: None
        """

        md_dict = copy.deepcopy(method_dict)

        # filter is set at construction
        md_dict['params']['filter'] = self.int_group

        # assign values to arguments we don't want to change
        if 'kT' in method_dict['th_params']:
            if method_dict['th_params']['kT'] == 'unchanged':
                md_dict['th_params']['kT'] = self.realtemp
                method_dict['th_params']['kT'] = self.realtemp
                
        if 'kT' in method_dict['params']:
            if method_dict['params']['kT'] == 'unchanged':
                md_dict['params']['kT'] = self.realtemp
                method_dict['params']['kT'] = self.realtemp

        if 'thermostat' in method_dict:
            therm_stat = getattr(hoomd.md.methods.thermostats, md_dict['thermostat'])(**md_dict['th_params'])
            self.method = getattr(hoomd.md.methods, md_dict['name'])(**md_dict['params'], thermostat=therm_stat)
        else:
            self.method = getattr(hoomd.md.methods, md_dict['name'])(**md_dict['params'])

        if self.integrator.methods:
            self.integrator.methods.pop()
            self.integrator.methods.append(self.method)

        # if not read from file, setup the new thermostat
        if not read_from_file:
            self.dict_thermostat_dof = method_dict

    def run_hoomd(self, sim_param, dt=None, log_partial_energies=False, additional_filters=None):
        """
        runs hoomd for as many time steps dumping the results

        :param sim_param: SimParameters object
        :param dt: time step
        :param log_partial_energies: whether to log the partial contributions to the potential energy.
        :param additional_filters: list of additional filters
        """

        # add drift removal default if not given
        if not self.drift_removal_dict:
            drift_trigger = hoomd.trigger.Periodic(int(sim_param.steps_sim/2))
            self.sim.operations.updaters.append(hoomd.md.update.ZeroMomentum(drift_trigger))

        # make sure that no forces are present
        self.integrator.forces.clear()
        # add all forces
        for pots in self.potentials.values():
            self.integrator.forces.append(pots)

        self.file_dict['write'] = self.file_write
        self.file_dict['restart'] = self.file_restart
        self.file_dict['simulation_data'] = self.file_sim
        self.file_dict['thermo_data'] = self.file_thermo
        self.file_dict['pickle'] = self.file_pickle

        # define logger 2 (adds thermodynamic quantities)
        # logger2 = hoomd.logging.Logger(categories=['scalar'])
        logger2 = hoomd.logging.Logger()
        logger2.add(self.sim, quantities=['timestep'])

        list_filters = [hoomd.filter.All()]
        if additional_filters is not None:
            for filts in additional_filters:
                list_filters.append(filts)
        for filts in list_filters:
            therm_prop = hoomd.md.compute.ThermodynamicQuantities(filter=filts)
            self.sim.operations.computes.append(therm_prop)
            for quant in sim_param.quantities_log:
                logger2.add(therm_prop, quantities=quant)
        # add partial energies if requested
        if log_partial_energies:
            for pots in self.potentials.values():
                logger2.add(pots, quantities=['energy'])

        # define logger 3 (adds simulation quantities)
        logger3 = hoomd.logging.Logger(categories=['scalar'])
        logger3.add(self.sim, quantities=['timestep', 'walltime', 'tps'])

        # define writers with triggers and loggers
        trig1 = hoomd.trigger.Periodic(sim_param.steps_write)
        lst_write = ['property', 'momentum']
        gsd_writer1 = hoomd.write.GSD(filename=self.file_write, trigger=trig1, mode='ab', dynamic=lst_write)
        trig2 = hoomd.trigger.Periodic(sim_param.steps_log)
        gsd_writer2 = hoomd.write.GSD(filename=self.file_thermo, trigger=trig2, mode='ab', filter=hoomd.filter.Null())
        trig3 = hoomd.trigger.Periodic(sim_param.steps_log)
        gsd_writer3 = hoomd.write.GSD(filename=self.file_sim, trigger=trig3, mode='ab', filter=hoomd.filter.Null())

        gsd_writer2.logger = logger2
        gsd_writer3.logger = logger3

        # add the standard writers to the simulations
        self.gsd_writers.append(gsd_writer1)
        self.gsd_writers.append(gsd_writer2)
        self.gsd_writers.append(gsd_writer3)

        if self.write_data:
            for wr in self.gsd_writers:
                self.sim.operations.writers.append(wr)

        # if dt is None use the default parameter set by units. Otherwise, change it
        if dt is not None:
            self.sim.operations.integrator.dt = dt

        # the restart_thermostat file will be eventually overwritten, so it is saved just in case
        if os.path.isfile(self.sysfile + '_restart_previous.gsd'):
            self.file_original_restart = self.sysfile + '_restart_thermostat_previous.json'
            self.file_dict['original_restart_thermostat'] = self.file_original_restart
            os.rename(self.sysfile + '_restart_thermostat.json', self.file_original_restart)

        # the _restart file will be eventually overwritten, so it is saved just in case
        self.file_original_restart = self.sysfile + '_restart_previous.gsd'
        self.file_dict['original_restart'] = self.file_original_restart
        os.rename(self.sysfile + '_restart.gsd', self.file_original_restart)

        # run the simulation
        self.sim.run(sim_param.steps_sim)

        # store thermostat degrees of freedom
        with open(self.sysfile + '_restart_thermostat.json', 'w') as tp:
            json.dump(self.dict_thermostat_dof, tp)

        # create the restart file
        hoomd.write.GSD.write(state=self.sim.state, filename=self.file_restart, mode='wb')
        # save a file encompassing all the files that are used in this simulation
        file_summary = self.sysfile + '_summary.txt'
        with open(file_summary, 'w') as fp:
            for key, value in self.file_dict.items():
                fp.write(key+' '+value+'\n')
            
    def minimize_fire(self, num_steps, steps_log, dict_params=None, def_filter=None):
        """
        runs hoomd.md.integrate.minimize_fire to relax the system

        :param num_steps: number of steps to iterate
        :param steps_log: number of steps at which quantities are logged
        :param dict_params: dictionary of parameters for FIRE (see HOOMD-Blue documentation)
        :param def_filter: group of particles to integrate
        """

        file_fire = self.sysfile + '_fire_restart.gsd'
        file_fire_thermo = self.sysfile + '_fire_therm.gsd'

        # if dict_params is None, then use the defaults
        if dict_params is None:
            dict_params = {}
            dict_params['dt'] = 0.001
            dict_params['force_tol'] = 1e-3
            dict_params['angmom_tol'] = 1e-3
            dict_params['energy_tol'] = 1e-8

        # use the default particles
        if def_filter is None:
            def_filter = self.int_group

        self.sim.operations.integrator = None
        self.integrator.forces.clear()
        self.integrator = hoomd.md.minimize.FIRE(**dict_params, integrate_rotational_dof=self.rot, rigid=self.rigid_body)
        self.integrator.methods.append(hoomd.md.methods.NVE(def_filter))
        for pots in self.potentials:
            self.integrator.forces.append(pots)
        self.sim.operations.integrator = self.integrator

        # log energies
        trig2 = hoomd.trigger.Periodic(steps_log)
        gsd_writer = hoomd.write.GSD(filename=file_fire_thermo, trigger=trig2, mode='ab', filter=hoomd.filter.Null())
        logger_fire = hoomd.logging.Logger(categories=['scalar'])
        logger_fire.add(self.sim, quantities=['timestep'])
        for pots in self.potentials:
            logger_fire.add(pots, quantities=['energy'])

        gsd_writer.logger = logger_fire
        self.sim.operations.writers.append(gsd_writer)

        while not self.integrator.converged:
            self.sim.run(num_steps)

        if self.integrator.converged:
            print('FIRE MINIMIZER HAS CONVERGED')

        hoomd.write.GSD.write(state=self.sim.state, filename=file_fire, mode='wb')
