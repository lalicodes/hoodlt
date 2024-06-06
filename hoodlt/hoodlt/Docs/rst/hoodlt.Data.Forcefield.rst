.. _ForcefieldInfo:

Forcefields
===========

Forcefields in HOODLT are excel files which are used to store the simulation parameters, charges, masses, types, etc.
that are needed to build and run simulations. See :mod:`hoodlt.Data.Forcefield.example-empty_forcefield.xlsx` for an
overview of what should at minimum exist in a forcefield. If you want to add more columns, for example, to
document where you got your forcefield parameters from, feel free to add that information in a new and unused column.
Each forcefield has 6 tabs: groups, nonbonded, bond, angle, dihedral, and improper. Let's break each of these down one
by one:

The Groups Tab
--------------

The groups tab in a forcefield has many columns: name, type, long name, molecular weight, and charge. It also has 3
columns mass, energy, and length which one should have one entry, but more about that later.

    name - this column should specify the types (except for rigid centers) in the system you wish to simulate

    type - this column serves as a reminder as to whether the corresponding type is united or a single atom

    long name - this column is simply a way to specify in more detail what the corresponding type really represents

    molecular weight  - this column should be the mass of the corresponding type in the mass unit you choose

    charge - this column should be the charge of the corresponding type in units of the fundamental charge

The next 3 columns specify your chosen units to express all numbers in the forcefield in. You should just have one
entry in each of these columns. A list of available unit systems is given below:

    mass - can only be amu

    length - can be angstrom or nanometer (nm)

    energy - can be ev or kj/mol

The Nonbonded Tab
-----------------

The nonbonded tab has many columns as well, and we break these down here:

    name - this column should specify the types (except for rigid centers) in the system you wish to simulate. If your system has nonbonded parameters which do not come from the mixing rules, add those pairs to this column and set the nonbonded parameters for the pair interaction explicitly in the other columns

    sigma - this column should have the nonbonded parameter :math:`\sigma`, expressed in the length unit chosen in the groups tab

    epsilon - this column should have the nonbonded parameter :math:`\epsilon`, expressed in the energy units chosen in the groups tab

    alpha - this column should have the nonbonded parameter :math:`\alpha`

    qualifier - in the cell below the one that says 'qualifier', write either 'lj' or 'force_shifted_lj' to specify the nonbonded potential you wish to use.

    combine sigma - in the cell below the one that says 'combine sigma' write either 'arithmetic' or 'geometric' to specify the type of mixing rules this forcefield will use for combining :math:`\sigma` parameters

    combine epsilon - in the cell below the one that says 'combine epsilon' write either 'arithmetic' or 'geometric' to specify the type of mixing rules this forcefield will use for combining :math:`\epsilon` parameters

    combine alpha - in the cell below the one that says 'combine alpha' write either 'arithmetic' or 'geometric' to specify the type of mixing rules this forcefield will use for combining :math:`\alpha` parameters

The Bond Tab
------------

Here is an explanation of all the bond tab's columns:

    name - type of the bond

    qualifier - bond potential that this type implments

    modify - reminder as to whether or not the bond parameters for this type should be able to change (for example: all CTR-CTR bonds should have modoify = 1)

    k - the bond constant for the harmonic potential, expressed in the appropriate energy unit per length unit squared

    r0 - the bond length at equilibrium for the harmonic potential, expressed in length units

    sigma - the bond parameter for a fene bond, expressed in the length unit

    epsilon - the bond parameter for a fene bond, expressed in the energy unit

keep in mind that if you have many different bond potentials in the same forcefield, each bond type must implement each potential

The Angle Tab
-------------

Here is an explanation of all of the tab's columns:

    name - type of the angle

    qualifier - angle potential that this type implments

    k - constant for the angle potential, expressed in energy units per square radian

    t0 - equlibrium angle, expressed in radians

The Dihedral Tab
----------------

Here is an explanation of the tab's columns

    name - type of the dihedral

    qualifier - dihedral potential that this type implements

    group1/2/3/4 - reminder of the 4 types that are a part of the dihedral potential

    k1/2/3/4 - the 4 constants for the dihedral potential, in energy unit

The Improper Tab
----------------

Here is an explanation of the tab's columns

    name - type of the improper

    qualifier - improper potential that this type implements

    group1/2/3/4 - reminder of the 4 types that are a part of the improper potential

    k - the constant for the improper potential, in energy unit per radian squared

    chi - equilibrium angle for the improper potential, in radians