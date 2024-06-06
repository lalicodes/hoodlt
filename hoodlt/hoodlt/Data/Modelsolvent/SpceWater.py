"""
:module: Water
:platform: Unix, Windows
:synopsis: Implements a class defining water solvents

.. moduleauthor:: Elizabeth Macias
"""


from __future__ import division
import numpy as np
from hoodlt.Data.Modelsolvent.SolventAbs import SolventAbs


class SpceWater(SolventAbs):
    """
    Defines SPC/E model H2O as a rigid body
    """

    def __init__(self, ff):
        """

        :param repeats: number of repeats on the chain, for example, decane should have repeats=10
        :param ff: the name of the forcefield used to build the solvent
        """

        super(SpceWater, self).__init__(1, ff, 4, 'SpceWater')


        # types (_Water, O, H)
        self.types = ['_H2O', 'OW', 'H']
        self.typeid = ['_H2O']*1 + ['OW']*1 + ['H']*2

        self.diameter = np.zeros(self.num_particles) + 0.08

        # masses
        self.mass[1] = self.ff_reader.get_molecular_weight('OW') # amu
        self.mass[2] = self.ff_reader.get_molecular_weight('H') # amu
        self.mass[3] = self.ff_reader.get_molecular_weight('H') # amu
        self.mass[0] = np.sum(self.mass[1:]) # amu


        # charges
        self.charge[1] = self.ff_reader.get_charge('OW')
        self.charge[2] = self.ff_reader.get_charge('H')
        self.charge[3] = self.ff_reader.get_charge('H')

        # defining constituent particle positions
        theta_spce = 109.47*np.pi/180 # H-O-H angle in radians
        roh = self.ff_reader.get_bond_r0('OW-H') # distance between H and O.
        pos_spce = [(0,0,0), (roh*np.sin(theta_spce/2), roh*np.cos(theta_spce/2),0),
            (roh*np.sin(-theta_spce/2), roh*np.cos(theta_spce/2), 0) ] # O, H, and H positions

        # center of mass of the rigid body
        mass_spce = np.array([self.mass[1],self.mass[2],self.mass[3]]) # O, H, H mass in amu
        rcm = np.sum(mass_spce*pos_spce,axis=0)/np.sum(mass_spce) # Center of mass of O,H,H combo.

        # particle positions
        self.position[0] = rcm -rcm # H2O
        self.position[1:] = pos_spce-rcm # O, H, and H positions relative to center of mass

        # rigid body stuff
        #self.moment_inertia[0] = np.diag(self.moment_of_inertia())
        temp_pos_spce = self.position[1:]
        I_spce = np.zeros((3,3),dtype=np.float64) # calc. moment of inertial, SUM m_i*(r^{vector}_i)**2
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if i == j:
                        I_spce[i,j] += mass_spce[k]*np.linalg.norm(temp_pos_spce[k])**2.0
                    I_spce[i,j] -= mass_spce[k]*temp_pos_spce[k,i]*temp_pos_spce[k,j]

        self.moment_inertia[0] = (I_spce[0,0],I_spce[1,1],I_spce[2,2])
        self.body = np.zeros(self.num_particles)

        #print("self.position", self.position)


    def get_vector(self):
        """
        See documentation in BasicSystemEntity
        """

        return self.position[1] - self.position[0]
