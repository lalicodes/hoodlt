"""
:module: AngAmuEvUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing units which are based off A, amu, and ev

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - This class is now THE units class for forcefields which use angstrom, amu, ev units
..                  - Added documentation
..                Alex Travesset <trvsst@ameslab.gov>, June 2021
..                  - Added angstrom to construction variable
..                Alex Travesset <trvsst@ameslab.gov>, March 2022
..                  - change the default time step
"""

from hoodlt.Data.Units.Boltzmann import Boltzmann as Bz
from hoodlt.Data.Units.PhysicalUnits import PhysicalUnits
from hoodlt.Data.Units.Permittivity import Permittivity as Pt
from hoodlt.Data.Units.PhysicalConstants import PhysicalConstants as pc


class AngAmuEvUnits(PhysicalUnits):
    """
    Unit System that uses Angstrom as the length unit, amu as the mass unit, and ev as the energy unit
    """

    def __init__(self):
        """
        Default temperature of the simulation is 386.817K
        """

        super(AngAmuEvUnits, self).__init__('AngAmuEv', 1, .071293, 1, 1/30,
                                            Bz.kb_ev, Pt.perm_ev_ang, .02,
                                            pc.ang_to_m, pc.amu_to_kg, pc.ev_to_joule)

        self.angstrom_to_construction = 1
