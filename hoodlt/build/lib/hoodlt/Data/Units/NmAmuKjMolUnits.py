"""
:module: NmAmuKjMolUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing units which are based off nanometer, amu, and kj/mol

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - Added angstrom to construction variable
..                  - changed dt
..                Elizabeth Macias <emacias@iastate.edu>, February 2022
..                  - changed the permittivity from kj_mol_ang to kj_mol_nm
"""

import numpy as np
from hoodlt.Data.Units.PhysicalUnits import PhysicalUnits
from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc
from hoodlt.Data.Units.Boltzmann import Boltzmann as Bz
from hoodlt.Data.Units.Permittivity import Permittivity as Pt


class NmAmuKjMolUnits(PhysicalUnits):
    """
    Unit System that uses nanometers as the length unit, amu as the mass unit, and kj/mol as the energy unit
    """

    def __init__(self):
        """
        Default temperature of the simulation is 386.817K
        """

        super(NmAmuKjMolUnits, self).__init__('NmAmuKjMol', 1, .071293, 1, 1/30,
                                               Bz.kb_kj_mol, Pt.perm_kj_mol_nm, .005/np.sqrt(pc.ev_to_kj_mol),
                                               pc.nm_to_m, pc.amu_to_kg, pc.kj_mol_to_joule)

        self.angstrom_to_construction = 0.1
