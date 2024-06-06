"""
:module: SolventAbs
:platform: Unix, Windows
:synopsis: Abstract class for solvents

..moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, July 2019
..                  - Removed unecessary functions
..               Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                 - moved align() to basic entity class
..                 - remove redundant fields
..                 - solvent_name is now set in basic entity class
"""

from hoodlt.Data.Modelconfigurations.BasicSystemEntity import BasicSystemEntity


class SolventAbs(BasicSystemEntity):
    """
    Abstract Class from which all solvents should inherit
    """

    def __init__(self, repeats, ff, num_particles, name):
        """
        :param repeats: number of repeat units in the solvent
        """

        super(SolventAbs, self).__init__(ff, num_particles, name)

        # number of repeats in the solvent, if necessary
        self.repeats = repeats

    def get_name(self):
        """
        Outputs the name of the solvent, in a format so that we can append it to the name of the file being dumped

        :return: a name (String) describing the solvent
        """

        name = self.name

        # append number of repeats if greater than 1
        if self.repeats > 1:
            name += ('-' + str(self.repeats))

        return name
