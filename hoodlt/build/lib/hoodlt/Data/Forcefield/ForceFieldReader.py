"""
:module: ForceFieldReader
:platform: Unix, Windows
:synopsis: Class to read values from the excel forcefields so that other HOODLT components have a simpler interface
to get forcefield parameters

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> April 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu>, June 2019
..                  - added impropers
..                  - updated documentation
..                  - updated the equals method
..                  - Nonbonded rcuts are now set pairwise, instead of using a system wide value
..                  - Different Mixing rules can now be applied to different nonbonded interaction parameters
..                  - Can override mixing rules for a given pair interaction by giving the pair explicitly in the
..                    nonbonded tab
..                  - Can scale nonbonded interaction energies by a factor by supplying the info in a dictionary
..                  - Added a TON of documentation
..                Jacob Austin <jaustin2@iastate.edu>, February 2020
..                  - Added functionality for special lj pairs and special coulomb pair interactions
..                Xun Zha <xzha@iastate.edu>, March 2021
..                  - Replaced pandas.ExcelFile(path).parse(sheet) by pandas.read_excel(path, sheet)
..                  - Removed __getstate__(), __setstate__() and __deepcopy__(), because ForceFieldReader does not
..                  - has attribute pandas.DataFrame anymore
..                Jianshe Xia <xiajs6075@iccas.ac.cn> July 2021
..                  - added the table force field for bond, angle and dihedral
..                Alex Travesset <trvsst@ameslab.gov> May 2022
..                  - simplified the file to make it compatible
"""

import importlib_resources
import pandas as pd
import numpy as np
class ForceFieldReader:
    """
    Reads parameters directly from a force field and provides a simple interface for other HOODLT components to get
    force field values.
    """

    def __init__(self, ff_name):
        """
        Creates the object

        :param ff_name: name of the excel file to read from, without the '_forcefield.xlsx'
        """

        # make an excelfile object
        self.name = ff_name
        dname = 'Data/Forcefield/' + self.name + '_forcefield.xlsx'
        ref = importlib_resources.files('hoodlt') / dname
        with importlib_resources.as_file(ref) as path:
            self.ff_path = path

    def get_list_attributes(self):
        """
        Get a list of all the attributes
        """
        xp = pd.ExcelFile(self.ff_path)
        return list(xp.sheet_names)

    def get_list_attributes_without_cutoff(self):
        """
        Get a list of all the attributes that do not require cutoffs
        """

        l_attr = self.get_list_attributes()
        l_attr.remove('groups')
        xp = pd.ExcelFile(self.ff_path)
        l_no_cutoff = []
        for attr in l_attr:
            xl = xp.parse(attr, keep_default_na=False)
            if not xl['has cutoff'][0]:
                l_no_cutoff.append(attr)

        return l_no_cutoff

    def get_potentials(self, attr):
        """
        Return the sheet attr as a list

        :param attr: attribute
        :return: list
        """
        xp = pd.ExcelFile(self.ff_path)
        xl = xp.parse(attr, keep_default_na=False)
        lst = xl['qualifier'].unique().tolist()
        # do not include potentials whose qualifier are classified as rigid
        if 'rigid' in lst:
            lst.remove('rigid')
        return lst

    def get_potentials_dict(self, attr, pot):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param attr: attribute
        :param pot: potential name
        return: dictionary with
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse(attr, keep_default_na=False)
        indh = self._search_df(df, 'name', pot)
        dh = df[indh]
        num_param = dh['numparameters'].item()

        dict_val = {}
        for int in range(1, num_param+1):
            key = dh['name'+str(int)].item()
            value = dh['param'+str(int)].item()
            dict_val[key] = value

        return dict_val

    def has_rcut(self, attr, pot):
        """
        Return true or false if the potential has cutoff

        :param attr: attribute
        :param pot: potential name
        return: bool
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse(attr, keep_default_na=False)
        indh = self._search_df(df, 'name', pot)
        dh = df[indh]
        var = bool(dh['has cutoff'].item())

        return var

    def get_potentials_params(self, attr, pot, int_name):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param attr: attribute
        :param pot: potential name
        :param int_name: interaction name
        return: dictionary
        """

        dict_val = self.get_potentials_dict(attr, pot)
        dict_params = {}
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse(attr, keep_default_na=False)

        indh = self._search_df(df, 'name', int_name)
        # this is necessary in case that names are repeated with different potentials
        dh = df[indh][df[indh]['qualifier'].str.contains(pot)]
        for ind, key in enumerate(dict_val):
            dict_params[key] = dh['param'+str(ind+1)].item()

        return dict_params

    def get_non_bonded(self, pot, intnm1, intnm2, alpha=1):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param pot: potential name
        :param intnm1: interaction 1
        :param intnm2: interaction 2
        :param alpha: rescaling factor
        return: dictionary
        """

        dict_dim = self.get_potentials_dict('nonbonded', pot)
        dict_params = {}
        # if it is a rigid body id there is no interaction
        if intnm1[0] == '_' or intnm2[0] == '_':
            for key in dict_dim:
                dict_params[key] = 0
        # then use mixing rules
        else:
            dict1 = self.get_potentials_params('nonbonded', pot, intnm1)
            dict2 = self.get_potentials_params('nonbonded', pot, intnm2)

            dict_rules = self._get_mixing_rules(pot)
            for key in dict1:
                val1 = dict1[key]
                val2 = dict2[key]

                rule = dict_rules[key]
                if rule == 'explicit':
                    dict_params[key] = self.get_potentials_params('nonbonded', pot, intnm1+'-'+intnm2)
                else:
                    dict_params[key] = self._apply_mixing_rules(val1, val2, rule)

                if dict_dim[key] == 'energy':
                    dict_params[key] = alpha*dict_params[key]

        return dict_params

    def get_non_bonded_rcut(self, pot, intnm1, intnm2, rcut):
        """
        Return a dictionary with keys parameters of the potential, values dimensions of the parameters

        :param pot: potential name
        :param intnm1: interaction 1
        :param intnm2: interaction 2
        :param rcut: cutoff in dimensionless units
        return: cut off in real units
        """

        dictv = self.get_non_bonded(pot, intnm1, intnm2)
        return list(dictv.values())[0]*rcut

    def get_exclusions(self):
        """
        Get the exclusions for the force fields
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('groups')
        filt = df['exclusions'].notnull()
        return df['exclusions'][filt].to_list()

    def get_type(self, particle):
        """
        returns the type (united or atom) of the given particle

        :param particle: the name of the particle
        :return: the type, atom or united
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('groups')
        dh = self._search_df(df, 'name', particle)
        return df[dh]['qualifier'].item()

    def get_chemical_name(self, particle):
        """
        Gets the chemical name of the particle, as specified in the forcefield

        :param particle: the name of the particle
        :return: the chemical of the chosen particle, from the groups tab
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('groups')
        dh = self._search_df(df, 'name', particle)
        return df[dh]['long name'].item()

    def get_molecular_weight(self, particle):
        """
        Gets the molecular weight of the particle

        :param particle: the name of the particle
        :return: the molecular weight of the particle, in the same units as in the forcefield
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('groups')
        dh = self._search_df(df, 'name', particle)
        return df[dh]['molecular weight'].item()

    def get_charge(self, particle):
        """
        Gets the charge of the particle

        :param particle: the name of the particle
        :return: the charge of the particle, exactly as it appears in the groups tab of the forcefield
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('groups')
        dh = self._search_df(df, 'name', particle)

        return df[dh]['charge'].item()

    def get_angle_t0(self, int_name, pot='harmonic'):
        """
        Gets the equilibrium angle for the angle type. Assumes a harmonic potential

        :param int_name: the name of the angle type
        :param pot: potential
        :return: the equilibrium  angle for this angle type
        """
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('angle')
        # check that the value exists
        if not self._search_df(df, 'name', int_name).any():
            print(int_name, 'is not a valid angle interaction')

        dict_v = self.get_potentials_params('angle', pot, int_name)

        return dict_v['t0']

    def get_bond_r0(self, int_name, pot='harmonic'):
        """
        Gets the equilibrium bond distance for the given bond type

        :param int_name: the name of the bond
        :param pot: potential
        :return: the equilibrium bond distance
        """
        # check that the value exists
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('bond')
        if not self._search_df(df, 'name', int_name).any():
            print(int_name, 'is not a valid bond interaction')

        dict_v = self.get_potentials_params('bond', pot, int_name)
        return dict_v['r0']

    def get_nbnd_sigma_single_particle(self, particle):
        """
        Gets the nonbonded sigma value for a given particle in the forcefield

        :param particle: the name of the particle
        :return: the sigma for the given particle
        """
        # check that the value exists
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('nonbonded')
        if not self._search_df(df, 'name', particle).any():
            print(particle, 'is not a valid particle type')

        dict_v = self.get_potentials_params('nonbonded', 'lj', particle)
        return dict_v['sigma']

    def get_units(self):
        """
        Creates a units object based on values in the groups tab in the forcefield

        :return: an instance of a units class
        """

        dname = 'Data/Forcefield/' + self.name + '_forcefield.xlsx'
        ref = importlib_resources.files('hoodlt') / dname
        with importlib_resources.as_file(ref) as path:
            xp = pd.ExcelFile(path)

        df = xp.parse('groups')

        mass = df['mass'].iloc[0]
        length = df['length'].iloc[0]
        energy = df['energy'].iloc[0]

        units_class_name = self._determine_units_class(mass, length, energy)

        pk_name = 'hoodlt.Data.Units.' + units_class_name
        nano_import = __import__(pk_name, fromlist=[''])
        units_obj = getattr(nano_import, units_class_name)()

        return units_obj

    def _search_df(self, df, key, value):
        """
        find the rows where df[key]=value

        :param df: dataframe
        :param key: key
        :param value: value
        :return: dataframe
        """

        if key in df.columns.values.tolist():
            mat = (df[key] == value)
            if mat.any():
                return mat
        print(key, value)
        raise ValueError('key or value does not exist in forcefield')

    def _get_mixing_rules(self, pot):
        """
        Gets one of the types of mixing rules defined in the forcefield

        :param pot: potential
        :return: dictionary with mixing rules
        """

        dict3 = self.get_potentials_dict('nonbonded', pot)
        dict_rules = {}
        xp = pd.ExcelFile(self.ff_path)
        df = xp.parse('nonbonded')
        dh = df[self._search_df(df, 'name', pot)]
        for ind, key in enumerate(dict3):
            dict_rules[key] = dh['combine'+str(ind+1)].item()

        return dict_rules

    @staticmethod
    def _apply_mixing_rules(value_1, value_2, rules):
        """
        Applies the mixing rules to 2 values

        :param value_1: one of the values
        :param value_2: the other value
        :param rules: a string, either geometric or arithmetic
        :return: the result of applying the mixing rules to the two values
        """
        # potentials are defined explicitly
        if rules == 'explicit':
            return rules
        # use combination rules
        if rules == 'arithmetic':
            val_final = np.average([value_1, value_2])
        elif rules == 'geometric':
            val_final = np.sqrt(value_1 * value_2)
        else:
            raise ValueError("Mixing rules must be arithmetic or geometric")

        return val_final

    def get_all_width(self, tab):
        """
        Gets all the different qualifiers on a given tab

        :param tab: the name of the tab
        :return: all the qualifiers on the tab
        """

        bond_tab = pd.read_excel(self.ff_path, tab)
        return bond_tab['width'].unique().tolist()[0]

    def get_bond_tfile(self, bond_type, qualifier='table'):
        """
        Gets the bond table file for the given bond type

        :param bond_type: the name of the bond
        :param qualifier: the potential that the bond type implements
        :return: the table file name for the bond
        """

        bond_tab = pd.read_excel(self.ff_path, 'bond')
        rows_w_qual_and_name = (bond_tab['qualifier'] == qualifier) & (
            bond_tab['name'] == bond_type)

        dname = 'Data/Forcefield/' + bond_tab[rows_w_qual_and_name]['filename'].tolist()[0]
        ref = importlib_resources.files('hoodlt') / dname
        with importlib_resources.as_file(ref) as path:
            filename = path

        return filename

    def get_angle_tfile(self, angle_type, qualifier='table'):
        """
        Gets the angle table file for the given angle type

        :param angle_type: the name of the angle
        :param qualifier: the potential that the angle type implements
        :return: the table file name for the angle
        """

        angle_tab = pd.read_excel(self.ff_path, 'angle')
        rows_w_qual_and_name = (angle_tab['qualifier'] == qualifier) & (
            angle_tab['name'] == angle_type)

        dname = 'Data/Forcefield/' + angle_tab[rows_w_qual_and_name]['filename'].tolist()[0]
        ref = importlib_resources.files('hoodlt') / dname
        with importlib_resources.as_file(ref) as path:
            filename = path

        return filename

    def get_dihedral_tfile(self, dihedral_type, qualifier='table'):
        """
        Gets the dihedral table file for the given dihedral type

        :param dihedral_type: the name of the dihedral
        :param qualifier: the potential that the dihedral type implements
        :return: the table file name for the dihedral
        """

        dihedral_tab = pd.read_excel(self.ff_path, 'dihedral')
        rows_w_qual_and_name = (dihedral_tab['qualifier'] == qualifier) & (
            dihedral_tab['name'] == dihedral_type)

        dname = 'Data/Forcefield/' + dihedral_tab[rows_w_qual_and_name]['filename'].tolist()[0]
        ref = importlib_resources.files('hoodlt') / dname
        with importlib_resources.as_file(ref) as path:
            filename = path

        return filename

    def _determine_units_class(self, mass, length, energy):
        """
        Returns the string name of the units class that corresponds to the units chosen in the forcefield

        :param mass: the mass unit given in the forcefield
        :param length: the length unit given in the forcefield
        :param energy: the energy unit given in the forcefield
        :return: a string which is the name of the units class corresponding to the class chosen in the forcefield
        """

        mass_formatted = self._format_mass_string(mass)
        length_formatted = self._format_length_string(length)
        energy_formatted = self._format_energy_string(energy)

        if mass_formatted != 'Cg' and length_formatted != 'Cg' and energy_formatted != 'Cg':
            units_class_name = length_formatted + mass_formatted + energy_formatted + 'Units'
        elif mass_formatted == 'Cg' and length_formatted == 'Cg' and energy_formatted == 'Cg':
            units_class_name = 'CourseGrainedUnits'
        else:
            raise ValueError(
                "Do not mix coarse grained and physical unit systems")

        return units_class_name

    @staticmethod
    def _format_mass_string(mass_string):
        """
        Formats the mass string given in the forcefield to be a string which is part of the name of a units class

        :param mass_string: entry in the mass column in the groups tab of the forcefield
        :return: a string formatted so it is part of the name of a units class
        """

        ms = mass_string.title()

        # there is really only 1 mass unit (amu) that could be used for simulations
        if ms == 'Amu':
            return ms
        elif ms == 'Dimensionless':
            return 'Cg'
        else:
            raise ValueError(
                "Invalid mass unit input, valid inputs include: amu, dimensionless")

    @staticmethod
    def _format_length_string(length_string):
        """
        Formats the length string given in the forcefield to be a string which is part of the name of a units class

        :param length_string: entry in the length column in the groups tab of the forcefield
        :return: a string formatted so it is part of the name of a units class
        """

        ls = length_string.title()

        if ls == 'A' or ls == 'Angstrom':
            return 'Ang'
        elif ls == 'Nm':
            return 'Nm'
        elif ls == 'Dimensionless':
            return 'Cg'
        else:
            raise ValueError(
                "Invalid length unit input, valid inputs include: A, Angstrom, nm, dimensionless")

    @staticmethod
    def _format_energy_string(energy_string):
        """
        Formats the energy string given in the forcefield to be a string which is part of the name of a units class

        :param mass_string: entry in the energy column in the groups tab of the forcefield
        :return: a string formatted so it is part of the name of a units class
        """

        es = energy_string.title()

        if es == 'Ev':
            return 'Ev'
        elif es == 'Kj/Mol':
            return 'KjMol'
        elif es == 'Dimensionless':
            return 'Cg'
        else:
            raise ValueError(
                "Invalid energy unit input, valid inputs include: ev, kj/mol, dimensionless")

    def __eq__(self, other):
        """
        Checks if one forcefield is equal to another. We define a forcefield to be equal to another if they read from
        the same forcefield excel file, and use the same input dictionaries

        :param other: the other forcefield reader
        :return: whether the two forcefield readers are equal or not
        """

        return self.name == other.name and self.ff_path == other.ff_path and self.special_rcut == other.special_rcut \
            and self.rcut == other.rcut and self.scale_epsilon == other.scale_epsilon
