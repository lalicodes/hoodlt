"""
:module: NmAmuEvUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing units which are based off nanometer, amu, and ev

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - Added angstrom to construction variable
..                Elizabeth Macias <emacias@iastate.edu>, February 2022
..                  - changed the permittivity from ev_ang to ev_nm
"""

from hoodlt.Data.Units.Boltzmann import Boltzmann as Bz
from hoodlt.Data.Units.PhysicalUnits import PhysicalUnits
from hoodlt.Data.Units.Permittivity import Permittivity as Pt
from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc


class NmAmuEvUnits(PhysicalUnits):
    """
    Unit System that uses nanometers as the length unit, amu as the mass unit, and ev as the energy unit
    """

    def __init__(self):
        """
        Default temperature of the simulation is 386.817K
        """

        super(NmAmuEvUnits, self).__init__('NmAmuEv', 1, .071293, 1, 1/30,
                                            Bz.kb_ev, Pt.perm_ev_nm, .01,
                                            pc.nm_to_m, pc.amu_to_kg, pc.ev_to_joule)

        self.angstrom_to_construction = 0.1
