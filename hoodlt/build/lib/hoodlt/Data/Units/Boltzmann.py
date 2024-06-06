"""
:module: Boltzmann
:platform: Unix, Windows
:synopsis: Defines the Boltzmann Constant

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June 2018
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
                    - class now refers to the constants in a static way, so we don't need to make objects to access the
                    constants
                    - added boltzmann constant in kj/mol/K
"""

from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc


class Boltzmann:
    """
    Defines the Boltzmann constant in different units
    
    """

    kb_joule = 1.380648e-23
    kb_kj_mol = kb_joule*1e-3*pc.moles_to_atoms
    kb_ev = 8.61733e-5
    kb_cg = 1.0
