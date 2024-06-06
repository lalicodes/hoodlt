.. _SimulationExample3:

Example 3: Water simulation with spce
=====================================

This example shows how to run a simulation with water molecules. In this case we use the SPC/E model

This simulation is a continuation of the one in :ref:`NcsExample5`.

.. code-block:: python

    from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
    from hoodlt.HOOMD.SimParameters import SimParameters
    from hoodlt.Data.Units.NmAmuKjMolUnits import NmAmuKjMolUnits

    # the system we simulate is 8 water molecules
    name = 'pure_solvent_cspce_1000_sSpceWater_ffOpls-Aa_uNmAmuKjMol'

    # It is of critical importance to use the same forcefield the system was built in,
    # its name is right in the name of the gsd file
    ff = 'opls-aa'

    # number of time steps in each window
    # for this type of simulation, there will only be one window, so this is the number of time steps to run the simulation for
    steps_wind = 1000

    # number of time steps per log to the .log file
    steps_log = 100

    quant_log = ['kinetic_temperature', 'pressure', 'potential_energy', 'volume']

    # number of time steps per write to the '_dump' file
    steps_write = steps_log

    # MPI variables
    num_procs = 1

    # thermodynamic conditions
    temp = 300

    # object which contains simulation parameters
    sim_params = SimParameters(steps_wind, steps_log, quant_log, steps_write, num_procs)

    # simulation object
    list_electrostatic = [64, 4, 0.8, 0]
    sim = HoomdSimulation(sysfile=name, rcut=3, ff_name=ff, cluster=False, temp_in_kelvin=temp, list_electrostatic=list_electrostatic)
    # this method call runs the simulation
    sim.run_hoomd(sim_params)


The files in your directory are similar to the ones in :ref:`SimulationExample1`.




