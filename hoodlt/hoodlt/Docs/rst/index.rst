.. HOODLT documentation master file, created by
   sphinx-quickstart on Thu Feb  9 14:04:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HOODLT
======

Purpose
=======

.. only:: html

.. |Citing-HOODLT| image:: https://img.shields.io/badge/cite-hoodlt-blue.svg
    :target: citing.html

.. |Downloads| image:: https://img.shields.io/badge/downloads-hoodlt-brightgreen.svg?style=flat
    :target: https://bitbucket.org/trvsst/hoodlt/wiki/Home

|Citing-HOODLT| |Downloads|

HOODLT is a program initially developed to compute free energies of crystals using
dynamical lattice theory (DLT) and thermodynamic integration. Over time, however, it has
been evolving to become a tool particularly focused to predictions in nanoparticle systems.
It includes:

* Inventory of simple and binary lattices including calculation of packing fractions.

* Tools to calculate free energies of general crystals.

* General tools to analyze simulations and configurations of crystalline assemblies.

* MD Simulations of general Nanocrystals with `HOOMD-blue <https://glotzerlab.engin.umich.edu/hoomd-blue/>`_.

* Implementation of the Orbifold Topological Model (OTM) for Nanocrystals.

* Available experimental data in nanoparticle superlattices.


Resources
=========

- `BitBucket Repository <https://bitbucket.org/trvsst/hoodlt>`_:
  Source code and issue tracker.
  Scripts to validate that Hoodlt performs accurate simulations.

Related Tools
=============

- `freud <https://freud.readthedocs.io/>`_:
  Analyze HOOMD-blue simulation results with the **freud** Python library.

- `foyer <https://foyer.readthedocs.io/>`_:
  Perform atom-typing and define classical molecular modeling force fields with **foyer**.

- `ParmEd <https://parmed.github.io/ParmEd/html/index.html>`_:
  Handle molecular dynamics parameter and topology files with **ParmEd**.

- `ovito <https://www.ovito.org/>`_:
  Visualize and analyze molecular dynamics simulation data with **ovito**.

.. rubric::
    Contents

.. toctree::
    :maxdepth: 1
    :caption: Guides

    getting.started


.. toctree::
   :maxdepth: 1
   :caption: Common Tasks

   how.to
   core.modules


.. toctree::
   :maxdepth: 1
   :caption: PYTHON API

   citing.rst

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   hoodlt.Data
   hoodlt.HOOMD
   hoodlt.Analysis

.. toctree::
   :maxdepth: 1
   :caption: Notes

   hoodlt.common_build_structure
   hoodlt.common_config_organization
   hoodlt.common_lattice_types

