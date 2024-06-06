.. _NcsExample1:

Example 1: Single Core With Ligand Chains
=========================================

We use :func:`hoodlt.Data.Modelconfigurations.ConfigurationBuilder` as shown below

.. code-block:: python

    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand as Hydrocarbon

    # this is the name of the forcefield you have chosen to use to build your configuration
    forcefield = "ncs-in-solvent"

    # this is the core object
    core = TO201(forcefield)

    # this is the ligand, which will be grafted to the core
    lig = Hydrocarbon(repeats=11, ff=forcefield)

    # this is the ConfigurationBuilder object, which is used to add nanoparticles, solvents, bonds, etc. to the configuration
    builder = ConfigurationBuilder()

    # here we have added a single nanoparticle to the configuration, centered at [0, 0, 0]
    builder.add_nc(core, [lig]*core.graft_num, [0, 0, 0])

    # Here we give our configuration an alias, to help us remember what arrangement of nanoparticles are in the configuration
    # this is completely optional
    builder.set_alias("Single")

    # now that the confiuration is built, we get our configuration back
    conf = builder.get_configuration()

    # here we save the configuration to a gsd file (so we can simulate it) and to a pickle file (for analysis later)
    save_config(conf)



After running this script, you will find some new files in your directory, their name and purpose is described below:

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_restart.gsd
    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_bonds.json
    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_restart.pickle

The names of the files might at first look confusing, but they describe what is in the system you just created.
Lets break down each component:

.. code-block:: bash

    'Au201' = you have a core which models a core with 201 atoms (both Au and S (thiol) needed to attach ligands)
    'Hydrocarbon-n11' = The Hydrocarbon ligand, which has n=11 repeating units in its chain
    '_cSingle' = you named your configuration Single
    '_pNcs-In-Solvent' = you used the parameters from the ncs-in-solvent forcefield to make your system
    '_uAngAmuEv' = your system uses angstrom, amu, and ev as the length, mass, and energy unit, respectively
    '_restart' = Indicates that this contains a configuration that can be used to run simulations

Now, let's discuss what the files are actually for. The first file, named

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_restart.gsd

is the file that is needed to run the simulation of your system, it is a snapshot that will be used by HOOMD to
begin the simulation. To run a simulation, just run a script similar to those described in the many
examples in :ref:`HOODLTSimulations` by entering

.. code-block:: bash

    python name-of-your-script.py

The next file, named

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_bonds.json

is a json file that contains some information about the simulation as well as modifications to the
force field. An important function is to bookkeep bonds among nanoparticles, where bond strength and distance is
different from the one read from the force field. The documenting parameters that are modified
from the force field, may also be needed for subsequent analysis.

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv.pickle

Is a file which stores the necessary information to convert a hoomd snapshot into a
:mod:`FunctionalizedConfiguration` object, see :ref:`HOODLTConfigOrg` for a detailed description.
That is, it takes the snapshots encoded in the .gsd file (the first file) and turn it back into
a configuration object, so it can be used after the simulation ends for analysis, to create other new
configurations or even to manipulate features of nanoparticles, solvent or substrates.

Armed with this knowledge, you can proceed to simulate your system and analyze the results.


