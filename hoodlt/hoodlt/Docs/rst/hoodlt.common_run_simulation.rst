.. _HOODLTrunsim:

How to Run a Simulation with HOOMD-Blue and Analyze the Results
===============================================================

HOODLT automatizes the process of running all atom and coarse-grained simulations involving all type of nanoparticles,
thus allowing the end-researcher to focus her effort on the science rather than on the technicalities of the simulation,
drastically reducing the possibility of unintended errors and also ensuring full reproducibility for any independent
researcher.

A successful simulations with HOODLT consists of a pipeline consisting of three stages:

- *Building a system*.
- *Running the system* using HOOMD-Blue.
- *Analyzing the outcome* of a simulations.

Each of them is described in different parts of the documentation.

**Building a system** is described in :ref:`HOODLTData` and many examples are provided. Explanation of force fields is given
in :ref:`ForcefieldInfo` and of units in :ref:`Units`.

**Running the system** is described in :ref:`HOODLTSimulations` and many examples are provided. Explanation of HOOMD-Blue
is provided in the  `HOOMD-blue web page <http://glotzerlab.engin.umich.edu/hoomd-blue/>`_  where there are excellent
tutorials.

**Analyzing the outcome** of a simulation is described in :ref:`Analysis` where many examples are provided.
Because the output data is basically hoomd output, for example a
`hoomd gsd file <https://gsd.readthedocs.io/en/stable/python-module-gsd.hoomd.html>`_ many software packages already
exist, and our module contains analysis for which no other software is available.

The different examples provide a sufficient broad number of examples that enable virtually any relevant situation for
nanoparticle simulations.







