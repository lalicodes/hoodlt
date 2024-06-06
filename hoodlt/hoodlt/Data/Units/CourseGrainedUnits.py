"""
:module: CourseGrainedUnits
:platform: Unix, Windows
:synopsis: Implements simple class for representing a course grained unit system

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, June 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - Made the class concrete instead of abstract
..                Xun Zha <xzha@iastate.edu>, June 2021
..                  - Added angstrom to construction variable
"""

from hoodlt.Data.Units.UnitsClass import UnitsClass
from hoodlt.Data.Units.Boltzmann import Boltzmann as Bz
from hoodlt.Data.Units.Permittivity import Permittivity as Pt


class CourseGrainedUnits(UnitsClass):
    """
    Class which represents all course grained unit systems
    """
    def __init__(self):

        super(CourseGrainedUnits, self).__init__('CourseGrained', 1, 1, 1, 1, Bz.kb_cg, Pt.perm_cg, 1/1000)

        self.angstrom_to_construction = 1
