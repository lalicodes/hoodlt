.. _SimulationExample7:

Example 7: Changing Integrator
==============================

The standard HOOMD integrator used for simulations is the NVT integrator, but for certain simulations, it may be
more useful to use a different integrator. HOODLT just defines a dictionary, and passes it to HOOMD.

The folowing script can be run after :ref:`SimulationExample3`.

.. code-block:: python

    from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
    from hoodlt.HOOMD.SimParameters import SimParameters
    from hoodlt.Data.Units.NmAmuKjMolUnits import NmAmuKjMolUnits

    # the system we simulate is 8 water molecules
    name = 'pure_solvent_cspce_1000_sSpceWater_ffOpls-Aa_uNmAmuKjMol'

    # It is of critical importance to use the same forcefield the system was built in,
    # its name is right in the name of the gsd file
    ff = 'opls-aa'

    units = NmAmuKjMolUnits()

    # number of time steps in each window
    # for this type of simulation, there will only be one window, so this is the number of time steps to run the simulation for
    steps_wind = 1000
    #steps_wind = int(5E3)

    # number of time steps per log to the .log file
    steps_log = 100
    #steps_log = int(steps_wind/1E3)

    quant_log = ['kinetic_temperature', 'pressure', 'potential_energy', 'volume']

    # number of time steps per write to the '_dump' file
    steps_write = steps_log

    # MPI variables
    num_procs = 1

    # thermodynamic conditions
    temp = 300
    Pa_in_one_atm = 101325
    P = 1 #bar
    P = P*0.986923 # atm
    P = P*Pa_in_one_atm # SI units
    P = P/units.pressure_construction_to_si # kJ/mol/nm3

    print("hoodlt pressure conversion:", 1/units.pressure_construction_to_si)

    # object which contains simulation parameters
    sim_params = SimParameters(steps_wind, steps_log, quant_log, steps_write, num_procs)

    # simulation object
    list_electrostatic = [64, 4, 0.8, 0]
    sim = HoomdSimulation(sysfile=name, rcut=3, ff_name=ff, cluster=False, temp_in_kelvin=temp, list_electrostatic=list_electrostatic)

    # change integrator
    NPT_int_dict = {'name': 'NPT', 'params': {'kT': 300*units.sim_temp_default/units.kelvin_temp, 'S': P, 'filter': 'unchanged', 'tau': 0.01, 'tauS': 0.1, 'couple': "xyz"}}

    sim.change_integration_method(NPT_int_dict)

    print("switched integration method")
    # this method call runs the simulation
    sim.run_hoomd(sim_params, log_partial_energies=True)

In this case we switched to the NPT integrator.
The parameters for the new integrator are in the following dictionary

.. code-block:: python

    NPT_int_dict = {'name': 'NPT', 'params': {'kT': 300*units.sim_temp/units.kelvin_temp, 'S': P_sim, 'filter': 'unchanged', 'tau': 0.01, 'tauS': 0.1, 'couple': "xyz"}}

with the same parameters as specified in HOOMD for the given integrator:

For documentation on all the types of integrators available in HOOMD, see `here <https://hoomd-blue.readthedocs.io/en/latest/index.html>`_.

and string 'unchanged' specifies that the parameter is unchanged from the original NVT.


