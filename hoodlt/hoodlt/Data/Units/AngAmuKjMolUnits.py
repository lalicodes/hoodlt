"""
:module: AngAmuKjMolUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing units which are based off A, amu, and kj/mol

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - This class is now THE units class for forcefields which use angstrom, amu, kj/mol units
..                  - Added documentation
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - Added angstrom to construction variable
"""

import numpy as np
from hoodlt.Data.Units.PhysicalUnits import PhysicalUnits
from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc
from hoodlt.Data.Units.Boltzmann import Boltzmann as Bz
from hoodlt.Data.Units.Permittivity import Permittivity as Pt


class AngAmuKjMolUnits(PhysicalUnits):
    """
    Unit System that uses Angstrom as the length unit, amu as the mass unit, and kj/mol as the energy unit
    """

    def __init__(self):
        """
        Default temperature of the simulation is 386.817K
        """

        super(AngAmuKjMolUnits, self).__init__('AngAmuKjMol', 1, .071293, 1, 1/30,
                                               Bz.kb_kj_mol, Pt.perm_kj_mol_ang, .05/np.sqrt(pc.ev_to_kj_mol),
                                               pc.ang_to_m, pc.amu_to_kg, pc.kj_mol_to_joule)

        self.angstrom_to_construction = 1
