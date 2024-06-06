"""
:module: ConfigurationBuilder
:platform: Unix, Windows
:synopsis: Builds configurations of NCs

.. moduleauthor:: Tommy Waltmann <tomwalt@iastate.edu> May 2019
.. history:
..                Tommy Waltmann <tomwalt@iastate.edu> July 2019
..                  - added support for binary reinitialized lattices
..                  - solvent can be added to non-rectangular boxes now
..                  - removed the old way of adding solvent
..                  - NCs are now added at the given position, regardless of their position before adding
..                Tommy Waltmann <tomwalt@iastate.edu> June 2019
..                  - added support for substrates
..                  - added solvent addition that doesn't change box size
..                  - generalized addition of ncs and building of lattices for ncs with more than one type of ligand
..                  - added a way to add more to configuration you want to keep building
..                  - added a way to add reinitialized nanoparticles to a configuration
..                  - cleaned up solvent addition slightly
..                  - forcefield reader can now determine the units
..                  - Documentation
..                  - masses are now calculated correctly in the solvent addition method
..                Xun Zha <xzha@iastate.edu> June 2021
..                  - changed parameters in method build_lattice() and in build_lattice_from_reinit()
..                Alex Travesset <trvsst@ameslab.gov> March 2022
..                  - Made it compatible with hoomd v3
..                  - Eliminated dependence on hoomd classes
"""
import copy
import copy as cp
import numpy as np
import numpy.linalg as la
from hoodlt.Data.Modelconfigurations.FunctionalizedParticle import FunctionalizedParticle
from hoodlt.Data.Modelconfigurations.FunctionalizedConfiguration import FunctionalizedConfiguration
from hoodlt.Data.Modelconfigurations.SolventBuilder import SolventBuilder
from hoodlt.Data.Modelconfigurations.LatticeFunctionalizedConfiguration import LatticeFunctionalizedConfiguration


class ConfigurationBuilder:
    """
    Class which is used to build all configurations
    """

    def __init__(self, init_conf=None):
        """

        :param init_conf: FunctionalizedConfiguration to add more system components to
        """

        self.box_def = False  # has the box been defined the box yet?
        self.ff_def = False  # has the forcefield been defined yet?

        if init_conf is not None:
            self.conf = cp.deepcopy(init_conf)
            self._check_ff(init_conf.ff_reader)

            # if its a lattice, the box should not change
            if isinstance(self.conf, LatticeFunctionalizedConfiguration):
                self.box_def = True
        else:
            self.conf = FunctionalizedConfiguration([], [], -1)  # start with an empty configuration, null box

    def add_nc(self, core_obj, list_lig_objs, position, orientation=None):
        """
        Adds an NC to the configuration

        :param core_obj: the core on the NC to be added
        :param list_lig_objs: a list of ligand objects, with length core_obj.graft_num
        :param position: position of the NC to be added
        :param orientation: orientation of nc given as a quaternion
        :return: NC index of the NC which was just added
        """

        # core must be centered at the origin
        if not np.allclose(core_obj.position[0], np.zeros(3)):
            raise ValueError('center cores must be placed at the origin')

        # make the nc
        nc = self._make_nc(core_obj, list_lig_objs)

        # if the new nc is not present add it to the list of distinct particles
        nm = nc.get_name()
        if not any(nm in dict_vals for dict_vals in self.conf.distinct_particles):
            dict_v = {nm: copy.deepcopy(nc)}
            self.conf.distinct_particles.append(dict_v)

        # shift to the correct position
        if orientation is not None:
            nc.rotate(orientation)

        nc.shift(-1 * nc.core.position[0])
        nc.shift(position)

        # update the configuration
        self.conf.particles.append(nc)
        self.conf.bondable_entities.append(nc)

        return len(self.conf.bondable_entities) - 1  # NC index can be used for bonding to other NCs in add_bond()

    def add_substrate(self, substrate_obj, origin):
        """
        Adds a substrate to the configuration
        :param substrate_obj: a substrate object.
        :param origin: position of the rigid center of the substrate after it is added to the configuration
        :return: None
        """

        self._check_ff(substrate_obj.ff_reader)
        subs = cp.deepcopy(substrate_obj)

        # if the new nc is not present add it to the list of distinct particles
        subs_name = substrate_obj.get_name()
        if not any(subs_name in dict_vals for dict_vals in self.conf.distinct_particles):
            dict_v = {subs_name: copy.deepcopy(substrate_obj)}
            self.conf.distinct_substrates.append(dict_v)

        subs.shift(-1 * subs.position[0])
        subs.shift(origin)

        self.conf.substrates.append(subs)
        self.conf.bondable_entities.append(subs)

        return len(self.conf.bondable_entities) - 1

    def add_reinit_nc(self, nc, position, orientation=None):
        """
        Adds the nc to the configuration at the given position

        :param nc: FunctionalizedParticle object
        :param position: the position to place the center of the FunctionalizedParticle at
        :param orientation: nc orientation given as a quaternion
        :return: NC index of the NC which was just added
        """

        nc_copy = cp.deepcopy(nc)

        if not np.allclose(nc_copy.core.position[0], np.zeros(3)):
            print(nc_copy.core.position[0])
            raise ValueError('center cores must be placed at the origin')

        # check for forcefield compatibility
        self._check_ff(nc_copy.core.ff_reader)

        # if the new nc is not present add it to the list of distinct particles
        if not any(nc_copy.get_name() in dict_v for dict_v in self.conf.distinct_particles):
            nc_bf = copy.deepcopy(nc_copy)
            nc_bf.body_frame()
            dict_nc = {nc_copy.get_name(): nc_bf}
            self.conf.distinct_particles.append(dict_nc)

        # shift and orient nc so its centered at the new position
        if orientation is not None:
            nc_copy.rotate(orientation)

        nc_copy.shift(-1*nc_copy.core.position[0])
        nc_copy.shift(position)

        # update the configuration
        self.conf.particles.append(nc_copy)
        self.conf.bondable_entities.append(nc_copy)

        return len(self.conf.bondable_entities) - 1

    def add_bond(self, index1, index2, bond_type):
        """
        Adds a bond between two NCs of NC index index1 and index2

        :param index1: the NC index of the first NC in the bond
        :param index2: the NC index of the second NC in the bond
        :param bond_type: the type of the bond
        :return:
        """

        pair = [index1, index2]
        pair_equ = [index2, index1]

        # if the bond is already present do not add
        for lst in self.conf.bonds:
            if pair in lst or pair_equ in lst:
                print('only one type of bond can be added between ', pair, ' skipping')
                return

        # add the bond type if it is not already present
        if bond_type not in self.conf.bonds_types:
            self.conf.bonds_types.append(bond_type)
            self.conf.bonds.append([])
            self.conf.bonds_dist.append(None)

        # add the bond
        ind = self.conf.bonds_types.index(bond_type)
        self.conf.bonds[ind].append(pair)

    def add_solvent(self, list_solv_objs, grid_spacing, box=None, pos_list=None):
        """
        Adds the given list of solvent objects to the box

        :param list_solv_objs: list of solvent objects to add to the box
        :param grid_spacing: defines the grid where solvent will be placed
        Increment this value if the simulation will not start after the solvent is added.
        :param box: (in construction units) for configurations where the box is not set.
        :param pos_list: option position of the solvent
        :return: Actual density the solvent was added at, in simulation units
        """

        if len(self.conf.solvent) > 0:
            print('Adding solvent to a configuration with solvent already')

        # set the box for non-lattice configurations
        if not self.box_def:
            if box is None:
                raise ValueError("Must set the density of solvent to add, if it is not a lattice")
            self.set_box(box)

        # prepare the box for adding solvent
        place_solvent = SolventBuilder(self.conf.box, grid_spacing, self.conf)

        # add the solvents
        for ind, solv_obj in enumerate(list_solv_objs):
            # if the new solvent is not present add it to the list of distinct particles
            if not any(solv_obj.name in dict_s for dict_s in self.conf.distinct_solvent):
                dict_solvent = {solv_obj.get_name(): cp.deepcopy(solv_obj)}
                self.conf.distinct_solvent.append(dict_solvent)
            solv = cp.deepcopy(solv_obj)
            self._check_ff(solv.ff_reader)
            if pos_list is None:
                place_solvent.insert(solv)
            else:
                place_solvent.add(solv, pos_list[ind])
            self.conf.solvent.append(solv)

        print('done placing solvent')

    def set_box(self, box_dim):
        """
        Sets the box dimensions, if they are not already set

        :param box_dim: new box dimensions, see FunctionalizedConfiguration.create_box for allowed inputs
        :return:
        """

        if self.box_def:
            raise ValueError("Box size is already set, you can't change it after it has been set")

        self.conf.create_box(box_dim)

        self.box_def = True

    def build_lattice(self, list_core_objs, list_lig_objs, lat_obj, orientation=None):
        """
        Builds a general lattice. Configuration should be empty when this method is called

        :param list_core_objs: list of core objects
        :param list_lig_objs: list of ligand objects
        :param lat_obj: lattice object
        :param orientation: orientation quaternion
        :return: None
        """

        # check arguments
        if len(self.conf.bondable_entities) > 0 or len(self.conf.solvent) > 0:
            raise ValueError("Can only build lattices from otherwise empty configurations")
        if not isinstance(list_core_objs, list):
            raise ValueError('Parameter "list_core_objs" needs to be a list of core object(s)')
        if not isinstance(list_lig_objs, list):
            raise ValueError('Parameter "list_lig_objs" needs to be a list of ligand object(s)')
        if len(list_core_objs) != len(list_lig_objs):
            raise ValueError('Number of given core objects should match the number of given ligand objects')

        # create list of ncs
        list_ncs = []
        for core_obj, lig_obj in zip(list_core_objs, list_lig_objs):
            list_ncs.append(self._make_nc(core_obj, lig_obj))

        # set the configuration, and the forcefield
        self.conf = LatticeFunctionalizedConfiguration(lat_obj, list_ncs)
        if orientation is not None:
            if isinstance(object, dict):
                self.conf.rotate(orientation['quaternions'], indices=orientation['indices'])
            elif isinstance(object, tuple):
                self.conf.rotate([orientation])
            else:
                txt1 = 'orientations is either a dictionary with keys=indices, quaternions and values a list of'
                txt2 = 'particle indices and the corresponding tuple rotations or just a single 4-tuple'
                raise ValueError(txt1+txt2)
        self.conf.ff_reader = list_core_objs[0].ff_reader
        self.ff_def = True
        self.box_def = True

    def build_lattice_from_reinit(self, list_ncs, lat_obj):
        """
        Builds a single component or binary lattice from copies of input nanoparticles

        :param list_ncs: list of FunctionalizedParticle object(s)
        :param lat_obj: lattice object
        :return: None
        """

        # check arguments
        if len(self.conf.bondable_entities) > 0 or len(self.conf.solvent) > 0:
            raise ValueError("Can only build lattices from otherwise empty configurations")
        if not isinstance(list_ncs, list):
            raise ValueError('Parameter "list_ncs" needs to be a list of NC object(s)')

        self._check_ff(list_ncs[0].core.ff_reader)
        for nc in list_ncs:
            self._check_ff(nc.ff_reader)

        # set the configuration, and the forcefield
        self.conf = LatticeFunctionalizedConfiguration(lat_obj, list_ncs)
        self.conf.ff_reader = list_ncs[0].core.ff_reader
        self.ff_def = True
        self.box_def = True

    def set_alias(self, alias):
        """
        Sets the alias of the configuration

        :param alias: the "name" for the configuration you built
        :return:
        """
        self.conf.alias = alias

    def get_configuration(self):
        """
        gets the configuration you have built, after setting the box based on the NCs that were added.

        :return: the configuration you built
        """

        if not self.box_def:
            # make sure the origin is included
            pos = np.array([[0.0, 0.0, 0.0]])
            for part in self.conf.particles:
                pos = np.concatenate((pos, part.core.position))
                for lig in part.ligands:
                    pos = np.concatenate((pos, lig.position))

            # this is the maximum separation among particles
            max_dist = np.amax(la.norm(pos[np.newaxis, :, :] - pos[:, np.newaxis, :], axis=2))
            # define the box as slightly larger than twice the maximum distance between particles
            c_fac = 2.2
            box_dim = [c_fac*max_dist, c_fac*max_dist, c_fac*max_dist]
            self.set_box(box_dim)

        return self.conf

    def _check_ff(self, ff_reader):
        """
        Ensures that the given forcefield reader is the same as the forcefield reader currently used to build the
        configuration

        :param ff_reader: forcefield reader object
        :return: None
        """

        if self.ff_def:
            if self.conf.ff_reader.name != ff_reader.name:
                print(self.conf.ff_reader.name, ff_reader.name)
                raise ValueError("Cannot use multiple force fields to instantiate a configuration")
        else:
            self.conf.ff_reader = ff_reader
            self.ff_def = True

    def _make_nc(self, core_obj, list_lig_objs):
        """
        Makes a nanoparticle from a core and a ligand object. Each object will be deep copied before being used to build
        the nanoparticle

        :param core_obj: core object
        :param list_lig_objs: list of ligand objects
        :return: a FunctionalizedParticle object which has the given core and ligands
        """

        # make sure only 1 one force field is used, and deep copy the inputs
        self._check_ff(core_obj.ff_reader)
        core = cp.deepcopy(core_obj)

        list_ligs = []
        for lig in list_lig_objs:
            self._check_ff(lig.ff_reader)
            list_ligs.append(cp.deepcopy(lig))

        # return the FunctionalizedParticle
        return FunctionalizedParticle(core, list_ligs)
