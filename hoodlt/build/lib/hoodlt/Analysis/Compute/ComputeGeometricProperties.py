"""
:module: ComputeGeometricProperties
:Platform: Windows, Unix
:synopsis: Computes the bonds, angles, dihedrals of a simulation

.. moduleauthor:: Elizabeth Macias <emacias@iastate.edu>, March 2022

"""
import numpy as np
import numpy.linalg as la
import gsd.hoomd
import matplotlib.pyplot as plt
from hoodlt.Analysis.Compute.ComputeStatistics import StatComputeTimeSeries
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions

class ComputeGeometricProperties:
    """
    this class computes bonds, angles, dihedrals for any range of time frames
    """
    def __init__(self, file_name):
        """

        :param file_name: file to analyze the trajectory
        """

        self.all_bonds, self.all_angles, self.all_dihedrals = [], [], []
        self.file_name = file_name + '_write.gsd'
        self.sys = gsd.hoomd.open(self.file_name, mode='rb')

    def compute_properties(self, frames=None, list_prop=[]):
        """
        :param frames: list of frames to compute properties
        :param list_prop: list containing the name of the snap properties to compute
        :return individual: a list that contains [properties], a list of lists containing each propertie's attributes, units dictionary
        """

        print("Computing properties ...")

        self.list_prop = list_prop
        snap = self.sys.read_frame(0)
        self.names = [[getattr(snap, prop).types[i] for i in getattr(snap, prop).typeid] for prop in self.list_prop]
        self.units = {'bonds': 'nm', 'angles': 'radians', 'dihedrals': 'degrees'}

        for frame in frames:
            print("Frame: ", frame)
            bonds, angles, dihedrals = self.compute_properties_by_frame(frame)
            self.all_bonds.append(bonds)
            self.all_angles.append(angles)
            self.all_dihedrals.append(dihedrals)

        g_types = [self.all_bonds, self.all_angles, self.all_dihedrals]
        self.properties = [g for g in g_types if g]

        return

    def plot_properties(self):

        for p, prop in enumerate(self.list_prop):

            #store molecule data
            data = np.stack(self.properties[p])
            type_names = self.names[p]

            time = list(range(len(data)))
            for i in range(int(np.size(data) / len(data))):
                plt.plot(time, data[:,i])
                plt.title(type_names[i])
                plt.ylabel(str(prop)+" ("+str(self.units[prop]+")"))
                plt.xlabel("time step")
                plt.savefig(type_names[i]+'_'+str(i))
                plt.close()
        return

    def compute_averages(self, start=0, end=-1):

        for p, prop in enumerate(self.list_prop):

            #store molecule data
            data = np.stack(self.properties[p])
            type_names = self.names[p]

            thefile = open(str(prop)+'_'+str(self.units[prop])+'_.txt', 'w')
            string = '%s\t'*3+'%s\n'
            thefile.write(string % tuple([str(prop)] + ['mean', 'standard_error', 'standard_deviation']))
            for i in range(int(np.size(data) / len(data))):
                compute = StatComputeTimeSeries(data[:,i], start, end)
                mean, se, sd = compute.compute_average(), compute.compute_std_error(), compute.compute_std_deviation()
                string = '%s\t'+'%.10f\t'*2 + '%.10f\n'
                thefile.write(string % tuple([type_names[i]] + [mean, se, sd]))
        return

    def compute_properties_by_frame(self, frame):
        """
        :param frame: frame to compute attributes
        :return individual: 3 lists
        """

        self.snap = self.sys.read_frame(frame)
        self.BF = BoxFunctions(self.snap.configuration.box)
        position = self._adjust_positions(self.snap.particles)
        group_types = [getattr(self.snap, prop).group for prop in self.list_prop]

        bonds, angles, dihedrals = [], [], []
        for group in group_types:
            g_pos = [[position[i] for i in grp] for j, grp in enumerate(group)] # all types of prop, 4 positions per types
            if len(group[0]) == 2:
                bonds = [self._compute_bond(pos) for pos in g_pos]
            if len(group[0]) == 3:
                angles = [self._compute_angle(pos) for pos in g_pos]
            if len(group[0]) == 4:
                dihedrals = [self._compute_dihedral(pos) for pos in g_pos]

        return bonds, angles, dihedrals

    def _compute_bond(self, pos):
        """
        Computes the bond given the positions of the 2 particles that form a bond group

        :param pos: an array of 2 position arrays
        :return: scalar
        """

        return np.linalg.norm(pos[1] - pos[0])

    def _compute_angle(self, pos):
        """
        Computes the angle given the positions of the 3 particles that form an angle group

        :param pos: an array of 3 position arrays
        :return: scalar
        """

        #compute theta
        v1 = pos[0] - pos[1] # first bond vector
        v2 = pos[2] - pos[1] # sec bond vector

        first_arg = np.dot(v1, v2)
        sec_arg = np.linalg.norm(v1)*np.linalg.norm(v2)

        return np.arccos(first_arg / sec_arg)

    def _compute_dihedral(self, pos):
        """
        Computes the dihedral angle given the positions of the 4 particles that form a dihedral group

        :param pos: an array of 4 position arrays
        :return: scalar
        """

        #compute theta
        u1 = pos[1] - pos[0] # first bond vector
        u2 = pos[2] - pos[1] # sec bond vector
        u3 = pos[3] - pos[2] # third bond vector

        u1crossu2 = np.cross(u1,u2)
        u2crossu3 = np.cross(u2,u3)

        first_arg = np.dot(u2, np.cross(u1crossu2, u2crossu3))
        sec_arg = np.linalg.norm(u2)*np.dot(u1crossu2, u2crossu3)

        theta = np.arctan2(first_arg, sec_arg)
        if theta <0:
            theta = (2*np.pi) + theta

        return theta*(180/np.pi)


    def _adjust_positions(self, particles):
        """
        helper function. This CHANGES POSITIONS of solvent molecules which are across the
        box

        :return: None
        """

        return [self.BF.unwrap(pos, image) for pos, image in zip(particles.position, particles.image)]



