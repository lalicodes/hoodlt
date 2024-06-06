"""
:module: CollectForceFieldData
:Platform: Windows, Unix
:synopsis: Gets data from a forcefield to be used in simulation analysis and calculations

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu>, May 2019
.. history:
"""

from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader


class CollectForceFieldData:
    """
    This class gets data from a chosen forcefield and scales the values if the user chooses
    """
    def __init__(self, ff):
        """

        :param ff: name of the forcefield to get the data from, without the "_forcefield.xlsx" at the end
        """

        self.ff_reader = ForceFieldReader(ff)

    def get_bond_k(self, bond_name):
        """
        Get the harmonic bond constant

        :param bond_name: name of the bond to get the constant for
        :return: the bond constant k, for the bond name
        """

        return self.ff_reader.get_potentials_params('bond', 'harmonic', bond_name)['k']

    def get_bond_r0(self, bond_name):
        """
        Get the r0 of the bond with the name bond_name

        :param bond_name: the name of the bond type
        :return: the equlibrium distance of the bond named bond_name
        """

        return self.ff_reader.get_potentials_params('bond', 'harmonic', bond_name)['r0']

    def get_units(self):
        """
        Returns the units used by this force field

        :returns: units object
        """

        return self.ff_reader.get_units()
