.. _HOODLTExplainUnits:

How Units Work in HOODLT ?
==========================

HOODLT units defines units through :mod:`PhysicalUnits`, a derived class of :mod:`UnitsClass` that
uses :mod:`PhysicalConstants`, where several constants are defined.

Four classes, each one implementing an explict unit system, are defined:

================= ===========  =====================   =========
      Class         Length              Mass            Energy
================= ===========  =====================   =========
  AngAmuEvUnits    Angstrom     :math:`M_{CH2}` amu      eV
AngAmuKjMolUnits   Angstrom     :math:`M_{CH2}` amu      kJ/mol
  NmAmuEvUnits     nanometer    :math:`M_{CH2}` amu      eV
 NmAmuKjMolUnits   nanometer    :math:`M_{CH2}` amu      kJ/mol
================= ===========  =====================   =========

Here :math:`M_{CH2}=\frac{1}{.071293} = 14.0266`.

.. note:: Class names for units must follow the name convention: **LengthMassEnergy** as shown in the table

The choice of units is made by the forcefield, and it is not possible to change it anywhere else.


Building a system: units
^^^^^^^^^^^^^^^^^^^^^^^^

When building a system units we use **construction units**: those are the units of the system in the table
with all coefficients being one. For example for :mod:`NmAmuKjMolUnits`, construction units are given as
:math:`1 \mbox{nm},  1\mbox{amu}, 1 \mbox{kJ}/\mbox{mol}`.

All files saved by HOODLT are in **simulation units**. When a hoodlt configuration
(a :mod:`FunctionalizedConfiguration` object) is converted to a HOOMD snapshot and saved,
conversion from construction to simulation units is applied.


Running a system: units
^^^^^^^^^^^^^^^^^^^^^^^

Note that **simulation units** (as defined by the table above) are different from
**construction units**. For example, for :mod:`NmAmuKjMolUnits` the
simulation units are :math:`1 \mbox{nm}, 14.0266\mbox{amu}, 1 \mbox{kJ}/\mbox{mol}`.

The simulations that are run with HOOMD-Blue use one of the units in the table. The corresponding gsd file has
already been saved with the desired units. The forcefield parameters are also read consistent with the
units, as they use :mod:`ScaledForceFieldReader`, a derived class from :mod:`ForceFieldReader`
that converts the forcefield parameters from construction to simulation units.


Why separation between construction and simulation units?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reason is that construction should be readable units, and simulation should be units that
make the relevant quantities of order 1. Therefore, this separation provides flexibility for
the most possible general situations.

A note about cut-offs units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters used for cut-offs in hoodlt are dimensionless. If rcut is the parameter
inbtroduced for the cut-off, each non-bonded interactions will have a cut off with dimensions
given as

.. math::

    \Lambda_{with dimensions} = \sigma \times rcut

where :math:`\sigma` is the diameter in simulation units (entirely determined by the force field).

Code snippets
^^^^^^^^^^^^^

The following snippet shows conversion of the time step units in an actual simulation:

 .. code-block:: python

    from hoodlt.HOOMD.HoomdSimulation import HoomdSimulation

    # force field (the opls forcefield uses AngAmuEvUnits)
    ff = "ncs-in-solvent"

    name='Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv'

    # simulation object
    sim = HoomdSimulation(sysfile=name, rcut=5, ff_name=ff)

    # show a list of all the untis defined
    lst = sim.units.simulation_units_defined
    print(lst)

    # dt in simulation units
    dt = sim.units.dt
    print('dt in simulation units', dt)

    # dt in construction units: Angstrom*(amu/eV)^{1/2}
    dt_construction_units =  dt*sim.units.time_simulation_to_construction
    print('dt in consruction units', dt_construction_units)

    # dt in SI units (second)
    dt_si = dt_construction_units*sim.units.time_construction_to_si
    print('dt in si units', dt_si)

    # dt in SI units can be obtained also
    dt_si_also = dt*sim.units.time_simulation_to_si
    print('dt in si units, obtained directly from simulation units', dt_si_also)

    # dt in fs (femtoseconds)
    dt_fs = dt_si*1e15
    print('in femtoseconds', dt_fs)


For this program to run, it will need the files:

.. code-block:: bash

    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_restart.gsd
    Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_bonds.json

as obtained, for example, in :ref:`HOODLTData`

If the predefined time step is too large or small, it can be changed before running the simulation

 .. code-block:: python

    # change the time step
    user_defined_dt = 0.001
    sim.run_hoomd(sim_parameters, dt=user_defined_dt)
