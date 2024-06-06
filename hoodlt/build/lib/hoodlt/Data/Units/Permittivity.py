"""
:module: Permittivity
:platform: Unix, Windows
:synopsis: Defines the Vacuum Permittivity Constant in different units

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
                    - class now refers to the constants in a static way, so we don't need to make objects to access the
                    constants
                    - added permittivity in C^2/(kj/mol*A)
..                Elizabeth Macias <emacias@iastate.edu>, February 2022
..                  - added permittivity in C^2/(kj/mol*nm)
..                  - added permittivity in C^2/(ev*nm)
"""

from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc


class Permittivity:
    """
    Defines the permititivity of vacuum :math:`\\epsilon_0`
    
    """

    perm_ev_ang = .005526349406
    perm_ev_nm = perm_ev_ang * 10 # added C^2/(ev*nm) 10Ang/1nm
    perm_kj_mol_ang = perm_ev_ang / pc.ev_to_kj_mol
    perm_kj_mol_nm = perm_kj_mol_ang * 10 # added C^2/(kj/mol*nm) 10Ang/1nm
    perm_cg = 1.0
