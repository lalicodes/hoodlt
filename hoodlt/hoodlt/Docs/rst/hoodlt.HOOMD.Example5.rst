.. _SimulationExample5:

Example 5: Running a binary lattice
===================================


.. code-block:: python

    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds
    from hoodlt.HOOMD.SimParameters import SimParameters

    # the system we simulate is a planar configuration of 3 NCs surrounding a central nc
    # this is similar to the system we created in building example 3
    name = 'cAu1072S-Hydrocarbon-n11+cAu201S-Hydrocarbon-n11_cMgzn2_l221_a100_ffDry-Ncs_uAngAmuEv'

    ff = 'dry-ncs'

    steps_sim = 800
    steps_log = 200
    quant_log = ['kinetic_temperature','potential_energy','pressure']
    steps_write = 200
    steps_hist=100

    num_procs = 1

    sim_params = SimParameters(steps_sim, steps_log, quant_log, steps_write, steps_hist, num_procs)

    # simulation object, note that it is different than in example 1
    sim = SimulationWithBonds(sysfile=name, rcut=5, ff_name=ff, temp_in_kelvin=387)

    # here we set the coefficents for both types of CTR-CTR bonds in our system
    # rinit is the equilibrium distance, defaulting to the current separation
    # lamda is a strength scaling factor, defaulting to 1
    sim.initialize('l_00p2', lamda=[0.02, 0.02, 0.02], log_hist='all')

    # this method call runs the simulation
    sim.run_hoomd_hist(sim_params, log_partial_energies=True)


Here, we decided that the default parameters for the bond coefficient were too large and
decided to soften them to a value of 0.02.

Then, if we need to continue the simulation further, we use the following code:

.. code-block:: python

    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds
    from hoodlt.HOOMD.SimParameters import SimParameters

    # the system we simulate is a planar configuration of 3 NCs surrounding a central nc
    # this is similar to the system we created in building example 3
    name = 'cAu1072S-Hydrocarbon-n11+cAu201S-Hydrocarbon-n11_cMgzn2_l221_a100_ffDry-Ncs_uAngAmuEv_l_00p2'


    ff = 'dry-ncs'

    steps_sim = 800
    steps_log = 200
    quant_log = ['kinetic_temperature','potential_energy','pressure']
    steps_write = 200
    steps_hist=100

    num_procs = 1

    sim_params = SimParameters(steps_sim, steps_log, quant_log, steps_write, steps_hist, num_procs)

    # simulation object, note that it is different than in example 1
    sim = SimulationWithBonds(sysfile=name, rcut=5, ff_name=ff, temp_in_kelvin=387)

    # this method call runs the simulation
    sim.run_hoomd_hist(sim_params, log_partial_energies=True)
