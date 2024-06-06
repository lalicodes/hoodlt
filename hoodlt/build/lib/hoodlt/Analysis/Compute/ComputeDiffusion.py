"""
:module: ComputeSolventProperties
:Platform: Windows, Unix
:synopsis: Computes diffusion properties

.. moduleauthor:: Elizabeth Macias <emacias@iastate.edu>, March 2022

"""
import freud
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from hoodlt.Analysis.Compute.ComputeWithTrajectory import ComputeWithTrajectory
from hoodlt.Data.Modelconfigurations.BoxFunctions import BoxFunctions
import scipy.stats as st

class ComputeDiffusion(ComputeWithTrajectory):
    """
    this class uses trajectory objects to compute diffusion properties
    """
    def __init__(self, trajectory):
        """

        :param trajectory: a Trajectory object
        """

        super(ComputeDiffusion, self).__init__(trajectory)

        self._adjust_positions()

    def compute_diffusion_properties(self, units=None, dt=0, name='PEO', ligand=False):
        """
        Uses 
        """

        self.ligand = ligand
        self.time = (np.array(self.traj.timesteps) - self.traj.timesteps[0])*dt
        self.units = units
        self.name = name

        print("Computing diffusion ...")

        Rg_all_frames = []
        e2e_all_frames = []
        com_all_frames = []
        for frame in self.traj.frames:
            print("Frame: ", frame)
            Rg, e2e, com = self.compute_properties_by_frame(frame)
            print(Rg, e2e, com)
            Rg_all_frames.append(Rg)
            e2e_all_frames.append(e2e)
            com_all_frames.append(com)

        self.Rg = np.array(Rg_all_frames)
        self.e2e = np.array(e2e_all_frames)
        self.com_all_frames = com_all_frames

        return 

    def compute_properties_by_frame(self, frame):
        """
        Uses 
        """

        config = self.traj.configurations[self.traj.frames.index(frame)]

        if self.ligand:
            solvent = []
            for particle in config.particles:
                solvent += [solv for solv in particle.ligands if solv.name==self.name]
        else:
            solvent = [solv for solv in config.solvent if solv.name==self.name]

        Rg, end2end, com_list = [], [], []
        for solv in solvent:
            com = np.average(solv.position, axis=0)
            Rg += [np.sqrt(np.mean(np.linalg.norm(solv.position - com, axis=1)**2))]
            end2end += [np.linalg.norm(np.subtract(solv.position[-1],solv.position[0]))]
            com_list += [com]

        return np.mean(Rg), np.mean(end2end), com_list

    def compute_msd(self):
        """
        returns: an array for each molecule that contains com for all frames, lit is one
        molecule with an array containing all com for all frames, but x-y-z 
        are not group yet untill after this call
        """

        com_by_mol = [np.split(lit, len(self.traj.frames)) for lit in np.hstack(self.com_all_frames)]
        rdiff_by_mol = [ (mol - mol[0]) for mol in com_by_mol]
        rdiff_squared = np.array([np.linalg.norm(mol, axis=1) for mol in  rdiff_by_mol])**2
        self.MSD = np.average(rdiff_squared, axis=0)  # average of each frame

        return 

    def compute_freud_msd(self, mode='direct'):
        """
        returns: an array for each molecule that contains com for all frames, lit is one
        molecule with an array containing all com for all frames, but x-y-z 
        are not group yet untill after this call
        """

        box = freud.box.Box(Lx=3, Ly=3, Lz=3)
        freud_msd = freud.msd.MSD(box=box, mode=mode)
        self.MSD = freud_msd.compute(positions=self.com_all_frames).msd

        return

    def compute_window_msd(self, mode='direct'):
        """
        returns: an array for each molecule that contains com for all frames, lit is one
        molecule with an array containing all com for all frames, but x-y-z 
        are not group yet untill after this call
        """

        box = freud.box.Box(Lx=3, Ly=3, Lz=3)
        freud_msd = freud.msd.MSD(box=box, mode=mode)
        self.MSD = freud_msd.compute(positions=self.com_all_frames).msd

        return 


    def e2e_plot(self):
        """
        :return: an array with histogram end-to-end bins and an array with the probability distribution
        """

        '''
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        plt.plot(self.time, self.e2e, marker='o', ls='None', label = 'simulation')
        plt.axhline(y=4.315, color='black', label = 'maximum length')
        plt.title(self.name+" end-to-end")
        plt.ylabel(r"$ \overline{ \vert r-r0 \vert } \quad (nm) $")
        plt.xlabel("t (s)")
        plt.legend()
        plt.savefig('end_to_end')
        plt.show()
        '''

        stats = np.mean(self.e2e), np.std(self.e2e)/len(self.e2e), np.std(self.e2e)
        print("end-to-end, mean, se, sd:", [round(s, 3) for s in stats])

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        a, x, b = plt.hist(self.e2e, bins=20, histtype='step')
        plt.close()
        a_norm = np.array(a/len(self.e2e))
        a_exc = np.nonzero(a_norm)[0]
        plt.plot(x[1:][a_exc], a_norm[a_exc], '-o')
        plt.title("PEO end-to-end distance")
        plt.ylabel("probability")
        plt.xlabel("end-to-end")
        plt.savefig(self.name+'_end_to_end_hist')
        plt.show()

        return x[1:][a_exc], a_norm[a_exc]

    def Rg_plot(self):
        """
        :return: radius of gyration
        """

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        plt.plot(self.time, self.Rg, marker='o', ls='None', label='simulation')
        plt.title(self.name+" Rg")
        plt.ylabel("Rg (nm)")
        plt.xlabel("t (s)")
        plt.legend()
        plt.savefig('Rg')
        plt.show()

        stats = np.mean(self.Rg), np.std(self.Rg)/len(self.Rg), np.std(self.Rg)
        print("Rg, mean, se, sd:", [round(s, 3) for s in stats])

        return

    def msd_plot(self):
        """
        :return: difference squared
        """

        # convert nm^2 to m^2
        msd_m = self.MSD*(self.units.length_construction_to_si**2)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        plt.plot(self.time, msd_m, marker='o', ls='None', label = 'simulation')
        plt.title(self.name+" Diffusion")
        plt.ylabel(r"$<|r-r0|^2>$"+" ("+r"$m^2$"+")")
        plt.xlabel("t (s)")
        plt.legend()
        plt.savefig('diffusion')
        plt.show()

        return 


    def diffusion_constant(self, start=0):
        """
        :return: diffusion constant
        """

        # convert nm^2 to m^2
        msd_m = self.MSD*(self.units.length_construction_to_si**2)

        def msd(t, D, b):

            return 6*D*t + b

        popt, pcov= curve_fit(msd, self.time[-start:], msd_m[-start:])
        print(popt)

        ylim = max(msd_m)*.03
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        plt.plot(self.time, msd(self.time, *popt), 'r-',label='fit: D= '+ str("{:.2E}".format(popt[0]) ))
        plt.plot(self.time, msd_m, marker='o', ms=3, ls='None', label = 'simulation')
        plt.ylim(-ylim, max(msd_m)+ylim)
        plt.title(self.name+" Diffusion")
        plt.ylabel(r"$<|r-r0|^2>$"+" ("+r"$m^2$"+")")
        plt.xlabel("t (s)")
        plt.legend()
        plt.savefig('diffusion_fit')
        plt.show()

        return tuple(popt)

    def _adjust_positions(self):
        """
        helper function. Unwraps particle positions.

        :return: None
        """

        for config in self.traj.configurations:
            bf = BoxFunctions(config.box)
            for i, solv in enumerate(config.solvent):
                solv.position = bf.unwrap(solv.position, solv.image)
            for particle in config.particles:
                for solv in particle.ligands:
                    solv.position = bf.unwrap(solv.position, solv.image)
