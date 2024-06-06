.. _Units:

Units
=====

The units are implicitly used throughout the process of building, simulating, and analyzing systems. They work with the
forcefield to provide proper scaling of forcefield values into and out of dimensionless simulation units.

Unit systems are built from length, mass, and energy units chosen by a forcefield. Each system has appropriate
constants defined which provide proper scaling into dimensionless simulation units. See the forcefield documentation at
:ref:`ForcefieldInfo` to see which length, mass, and energy units are currently approved for use in a forcefield.

To add a unit system, simply create a class similar to :mod:`hoodlt.Data.Units.AmuAngEvUnits`, and alter
the method :func:`get_units` in :mod:`hoodlt.Data.Forcefield.ForceFieldReader` to correctly recognize the new unit
system. Keep in mind the new units class must set :math:`k_B, \epsilon_0` in the appropriate length and energy unit
for the new unit system.

Full list of unit systems currently available:

.. toctree::
   :maxdepth: 1

   hoodlt.Data.Units.list_all

Units classes can also be used separately from the forcefield to provide quick unit conversions from simulation numbers
to SI and other conventionally used units. Aside from that, there are many universal constants in HOODLT below

Universal Constants:

.. toctree::
   :maxdepth: 1

   hoodlt.Data.Units.Constants
