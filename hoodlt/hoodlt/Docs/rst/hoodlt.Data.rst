.. _HOODLTData:

Construction of Nanocrystal Systems
===================================

Introduction
------------

HOODLT is suitable for running simulations and proceeds according to these three steps:

- *Building a system*.
- *Running the system* using HOOMD-Blue.
- *Analyze the outcome* of the simulations.

The purpose of the modules here address the first part: **Building a system**. We need to:

- **Select NCs** with one or many cores, ligands, solvent and possibly a substrate.
- **Choose a force field** describing how the different atoms/particles constituting the system interact.
- **Define a particular structure** for the system: A cluster, a lattice, etc...

The purpose of these modules is to build the necessary objects to make a system and then write it to a .gsd file. This
file will be used to run MD simulations with HOOMD at a later stage.

Building all configurations is routed through one class: :mod:`hoodlt.Data.Modelconfigurations.ConfigurationBuilder`.
You will first begin by creating a :mod:`hoodlt.Data.Modelconfigurations.ConfigurationBuilder` object, which represents
an empty system. By calling the methods of this class, you will add components to your system, each defined using the
same forcefield. You can query the system you are building by calling the :func:`get_configuration` method.

The call to :func:`save_config` function defined in :mod:`hoodlt.Data.Modelconfigurations.Saver` module will save
the system into the .gsd file that runs a simulation using HOOMD.

All available cores are described in :ref:`ListCores`

All available ligands are described in :ref:`ListLigands`.

All available solvents are described in :ref:`ListSolvents`.

All available substrates are described in :ref:`ListSubstrates`

Therefore, you make core, ligand, solvent, and substrate objects which have been predefined in HOODLT, use the methods
of the :mod:`ConfigurationBuilder` class to make any configuration you wish, and then save your configuration so that
it can be simulated in HOOMD.

Here are some examples:

.. toctree::
   :maxdepth: 1

   hoodlt.Data.Example1
   hoodlt.Data.Example2
   hoodlt.Data.Example3
   hoodlt.Data.Example4
   hoodlt.Data.Example5
   hoodlt.Data.Example6
   hoodlt.Data.Example7
   hoodlt.Data.Example8
   hoodlt.Data.Example9


The packages in which the core, ligand, solvent, and substrate objects are defined are listed here:

.. toctree::
   :maxdepth: 1

   hoodlt.Data.Methods

If you would like to add your own core, ligand, solvent, or substrate to HOODLT, some conventions you should follow are located in :ref:`HOODLTConfigOrg`

Explanation of the Forcefield and Units is located here:

.. toctree::
   :maxdepth: 1

   hoodlt.Data.Forcefield
   hoodlt.Data.Units