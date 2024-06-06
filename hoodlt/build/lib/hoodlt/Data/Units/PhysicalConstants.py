"""
:module: PhysicalConstants
:platform: Unix, Windows
:synopsis: Class that holds many important physical constants

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
"""


class PhysicalConstants:
    """
    Class which stores many physical constants (unit conversions), the boltzmann and permittivity constants are defined
    in their respective classes
    """

    # fundamental conversions
    moles_to_atoms = 6.02214076e23
    pascal_to_atm = 9.86923e-6
    atm_to_bar = 1.01325
    m3_to_l = 1000
    ang_to_m = 1e-10
    nm_to_m = 1e-9
    ev_to_joule = 1.6021766208e-19  # also e to C
    amu_to_kg = 1.660539040e-27

    # derived conversions
    kj_mol_to_joule = 1e3/moles_to_atoms
    ev_to_kj_mol = ev_to_joule / kj_mol_to_joule
