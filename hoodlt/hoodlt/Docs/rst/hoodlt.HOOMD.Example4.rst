.. _SimulationExample4:

Example 4: Serial Potential of Mean Force
=========================================

The 2-body potential of mean force can tell us a lot about how a system behaves, and HOODLT has a simulation routine that
is tailored to facilitate the calculation of the PMF from its output data. The following simulation creates output data
that can be used by the analysis tools in HOODLT to easily calculate the potential of mean force via the weighted histogram
analysis method.

The Weighted Histogram Analysis Method `(WHAM) <http://www.alchemistry.org/wiki/Weighted_Histogram_Analysis_Method>`_

.. code-block:: python

    import numpy as np
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelconfigurations.Saver import load_config
    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelnanoparticles.AuS import AuS
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand
    from hoodlt.Data.ProcessConfigurations.Squeeze import Squeeze
    from hoodlt.Analysis.Collect.ReInitHelper import reinit_config
    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds
    from hoodlt.HOOMD.SimParameters import SimParameters

    forcefield = 'opls_dry-ncs_mix'

    name = 'cAu4033S-Hydrocarbon-n17_cSqueeze_pOpls_Dry-Ncs_Mix_uAngAmuEv_restart'

    conf_nc = reinit_config(name, name)

    sep_list = np.arange(70,111,1,dtype=int)

    for separation in sep_list:
    # Building configuration
        builder_pair = ConfigurationBuilder()
        half_separation = separation/2
        print("Separation of NC centers is", separation)

        ind1 = builder_pair.add_reinit_nc(conf_nc.particles[0], [-half_separation, 0, 0])
        ind2 = builder_pair.add_reinit_nc(conf_nc.particles[0], [half_separation, 0, 0])

        print('setting bond')
        builder_pair.add_bond(ind1, ind2, 'CTR-CTR1')

        l_box = float(2.5*separation)
        print('setting box size as', l_box)
        builder_pair.set_box(l_box)

        pair_alias = 'SphereSphere_'+str(separation)

        builder_pair.set_alias(pair_alias)

        conf = builder_pair.get_configuration()
        save_config(conf)

        # simulation parameters
        steps_sim = 6000
        steps_log = 500
        quant_log = ['kinetic_temperature','potential_energy','pressure']
        steps_write = 2000
        steps_hist= 1000
        num_procs=1
        sim_params = SimParameters(steps_sim, steps_log, quant_log, steps_write, steps_hist, num_procs)

        name_file = 'cAu4033S-Hydrocarbon-n17_c'+pair_alias+'_pOpls_Dry-Ncs_Mix_uAngAmuEv'

        # simulation object
        sim = SimulationWithBonds(sysfile=name_file, rcut=5, ff_name=forcefield, cluster=True, temp_in_kelvin=387)

        # here we set the coefficents for both types of CTR-CTR bonds in our system
        # rinit is the equilibrium distance, defaulting to the current separation
        # lamda is a strength scaling factor, defaulting to 1
        sim.initialize('l_00p2', rinit=[float(separation)], lamda=[0.02], log_hist='average')

        # this method call runs the simulation
        sim.run_hoomd_hist(sim_params, log_partial_energies=True)

This simulation as stated here will take significant amount of time. It is advised that it is run in a gpu cluster. This
is what the flag cluster=True implies.

The simulations can be continued for as many times and as many time steps as necessary. For example:

.. code-block:: python

    import numpy as np
    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds
    from hoodlt.HOOMD.SimParameters import SimParameters

    # forcefield
    forcefield = 'opls_dry-ncs_mix'

    # simulation parameters
    steps_sim = 2_000_000
    steps_log = 10_000
    quant_log = ['kinetic_temperature','potential_energy','pressure']
    steps_write = 500_000
    steps_hist= 10_000
    num_procs=1
    sim_params = SimParameters(steps_sim, steps_log, quant_log, steps_write, steps_hist, num_procs)

    sep_list = np.arange(70,111,1,dtype=int)

    for separation in sep_list:

        pair_alias = 'SphereSphere_'+str(separation)
        name_file = 'cAu4033S-Hydrocarbon-n17_c'+pair_alias+'_pOpls_Dry-Ncs_Mix_uAngAmuEv'+'_l_00p2'

        # simulation object
        sim = SimulationWithBonds(sysfile=name_file, rcut=5, ff_name=forcefield, cluster=True, temp_in_kelvin=387)

        sim.run_hoomd_hist(sim_params,log_partial_energies=True)


In addition to all the files discussed in :ref:`SimulationExample1` and :ref:`SimulationExample3`, there will again be
some more files in your directory after the simulation finishes, each containing output data

.. code-block:: bash

    cAu4033S-Hydrocarbon-n17_cSphereSphere_XX_pOpls_Dry-Ncs_Mix_uAngAmuEv_l_00p2_hist.gsd

where XX is each one of the distances. These files are necessary to generate the pmf.
