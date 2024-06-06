"""
:module: ScaledForceFieldReader
:platform: Unix, Windows
:synopsis: Class to read values from the excel forcefields, and then scale the values by the units system,
so that those numbers can be immediately inputted into HOOMD

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> April 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - units are now taken from ff instead of as argument
..                  - added method for implementing impropers
..                  - Added all the documentation
..                  - Nonbonded rcuts are now set pairwise, instead of using a system wide value
..                  - Different Mixing rules can now be applied to different nonbonded interaction parameters
..                  - Can override mixing rules for a given pair interaction by giving the pair explicitly in the
                      nonbonded tab
..                  - Can scale nonbonded interaction energies by a factor by supplying the info in a dictionary
..                Jacob Austin <jaustin2@iastate.edu>, February 2020
..                  - Added functionality for special lj pairs and special coulomb pair interactions
..                Jianshe Xia <xiajs6075@iccas.ac.cn> July 2021
..                  - added the table force field for bond, angle and dihedral
..                Alex Travesset <trvsst@ameslab.gov>
..                  - fixed a bug that would have given an error if sim/cons units used different length units
"""

from hoodlt.Data.Forcefield.ForceFieldReader import ForceFieldReader
class ScaledForceFieldReader(ForceFieldReader):
    """
    Reads parameters directly from a force field and provides a simple interface for other HOODLT components to get
    force field values. All returned values are scaled by the chosen unit system
    """

    def __init__(self, ff_name):
        """

        :param ff_name: name of the forcefield to read from
        """

        super(ScaledForceFieldReader, self).__init__(ff_name)

        self.units = self.get_units()

    def get_potentials_params(self, attr, pot, int_name):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param attr: attribute
        :param pot: potential name
        :param int_name: interaction name
        return: dictionary
        """

        dict_v = super(ScaledForceFieldReader, self).get_potentials_params(attr, pot, int_name)
        dict_dim = super(ScaledForceFieldReader, self).get_potentials_dict(attr, pot)
        new_dict = {}
        for key in dict_v:
            uname = dict_dim[key]
            fac = getattr(self.units, uname+'_construction_to_simulation')
            new_dict[key] = dict_v[key]*fac

        return new_dict

    def get_non_bonded(self, pot, intnm1, intnm2, alpha=1):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param pot: potential name
        :param intnm1: interaction 1
        :param intnm2: interaction 2
        :param alpha: rescaling factor
        return: dictionary
        """
        dict_v = super(ScaledForceFieldReader, self).get_non_bonded(pot, intnm1, intnm2, alpha=alpha)
        dict_dim = super(ScaledForceFieldReader, self).get_potentials_dict('nonbonded', pot)
        new_dict = {}
        for key in dict_v:
            uname = dict_dim[key]
            fac = getattr(self.units, uname + '_construction_to_simulation')
            new_dict[key] = dict_v[key] * fac

        return new_dict

    def get_non_bonded_rcut(self, pot, intnm1, intnm2, rcut):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param pot: potential name
        :param intnm1: interaction 1
        :param intnm2: interaction 2
        :param rcut: cutoff in dimensionless units
        return: cut off in real units
        """

        c_fact = self.units.length_construction_to_simulation
        val = super(ScaledForceFieldReader, self).get_non_bonded_rcut(pot, intnm1, intnm2, rcut)

        return val*c_fact

    def get_molecular_weight(self, particle):
        """
        Gets the molecular weight of the particle
        :param particle: the name of the particle
        :return: the molecular weight of the particle in dimensionless units
        """

        m_fact = self.units.mass_construction_to_simulation
        return super(ScaledForceFieldReader, self).get_molecular_weight(particle) * m_fact

    def get_charge(self, particle):
        """
        Gets the charge of the particle
        :param particle: the name of the particle
        :return: the charge of the particle in dimensionless units
        """

        c_fact = self.units.charge_construction_to_simulation
        return super(ScaledForceFieldReader, self).get_charge(particle) * c_fact

    def get_bond_r0(self, int_name, pot='harmonic'):
        """
        Gets the equilibrium bond distance for the given bond type

        :param int_name: the name of the bond
        :param pot: potential
        :return: the equilibrium bond distance
        """
        c_fact = self.units.length_construction_to_simulation
        return super(ScaledForceFieldReader, self).get_bond_r0(int_name, pot) * c_fact

    def get_nbnd_sigma_single_particle(self, particle):
        """
        Gets the nonbonded sigma value for a given particle in the forcefield

        :param particle: the name of the particle
        :return: the sigma for the given particle
        """
        c_fact =  self.units.length_construction_to_simulation
        return super(ScaledForceFieldReader, self).get_nbnd_sigma_single_particle(particle) * c_fact
