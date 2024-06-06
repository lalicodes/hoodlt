.. _SimulationExample6:

Example 6: Energy minimization with FIRE
========================================

Oftentimes it is convenient to obtain the minimum of energy. HOOMD provides FIRE
see `here <https://hoomd-blue.readthedocs.io/en/latest/module-md-minimize.html?highlight=FIRE#hoomd.md.minimize.FIRE>`_.

which can be directly called from HOODLT. Here is the example script:

.. code-block:: python

    from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation
    from hoodlt.HOOMD.SimParameters import SimParameters

    # the system we simulate is a single Au140S Hydrocarbon 8 repeat nanoparticle
    # this is very similar to the system we created in building example 1
    name = 'Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv'

    # It is of critical importance to use the same forcefield the system was built in,
    # its name is right in the name of the gsd file
    ff = "ncs-in-solvent"

    # number of time steps in each window
    # for this type of simulation, there will only be one window, so this is the number of time steps to run the simulation for
    num_steps = 1000

    # number of time steps per log to the .log file
    steps_log = 50

    # simulation object
    sim = HoomdSimulation(sysfile=name, rcut=5, ff_name=ff, temp_in_kelvin=387)

    # this method call runs the simulation
    sim.minimize_fire(num_steps, steps_log)


