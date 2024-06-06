.. _SimulationExample1:

Example 1: Basic Simulation
===========================

The most basic type of simulation that can be run is just simulating a system without CTR-CTR bonds for a given number
of timesteps, as shown below:

.. code-block:: python

    from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
    from hoodlt.HOOMD.SimParameters import SimParameters

    # temperature in Kelvin
    temp = 387

    # the system we simulate is a single Au140S Hydrocarbon 8 repeat nanoparticle
    # this is the system we created in building example 1
    name = 'Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv'

    # It is of critical importance to use the same forcefield the system was built in,
    # its name is right in the name of the gsd file
    ff = "ncs-in-solvent"

    # number of time steps in each window
    # for this type of simulation, there will only be one window, so this is the number of time steps to run the simulation for
    #steps_sim = 4000000
    steps_sim = 400

    steps_log = 200
    quant_log = ['kinetic_temperature','potential_energy','pressure']
    steps_write = 200
    steps_hist = 100

    # MPI variables
    num_procs = 1

    # object which contains simulation parameters
    sim_params = SimParameters(steps_sim, steps_log, quant_log, steps_write, num_procs)

    # simulation object
    sim = HoomdSimulation(sysfile=name, rcut=5, ff_name=ff, temp_in_kelvin=temp)

    # this method call runs the simulation
    sim.run_hoomd(sim_params, dt=0.02, log_partial_energies=True)



Before you ran this script, you should have built the configuration, see :ref:`HOODLTData` for examples, and these files
should have been in your directory:

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_restart.gsd
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_bonds.json
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_restart.pickle

After you run this script, there will be some more files in your directory, each containing output data

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_summary.txt
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_restart_previous.gsd
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_write.gsd
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_therm.gsd
    Au201-Hydrocarbon-n11_cSingle_ffNcs-In-Solvent_uAngAmuEv_sim.gsd

The '.txt' file contains a list of all the files associated with the run.

The '_restart_previous.gsd' is the gsd restart file used at the beginning of the simulation.
It is saved in case you want to start the simulation from the very beginning. The new
'_restart.gsd' gsd file contains the state of the system at the end of the simulation.

The '_write.gsd' gsd file contains snapshots of the simulation at intermediate times.

The '_therm.gsd' gsd files contains the actual quantities measured during the simulation

The '_sim.gsd' gsd files contains information about the simulation itself, such as timesteps per second,
duration of the simulation.
