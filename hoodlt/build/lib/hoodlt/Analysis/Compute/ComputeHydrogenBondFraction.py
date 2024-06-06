"""
:module: ComputeHydrogenBonds
:Platform: Windows, Unix
:synopsis: Computes Hydrogen Bonds fracctions by solving system of equations

.. moduleauthor:: Elizabeth Macias <emacias@iastate.edu>, November 2022
.. history:
..                
..                Elizabeth Macias <emacias@iastate.edu>, November 2022
..                  - adapted code from Aqpolypy (Author: Alex Travesset <trvsst@ameslab.gov>)
..                  - made it compatible with hoomd v3
..                  - added capability to fit for change-in-energy and change-in-entropy parameters
..                  - added capability to plot numerical result
..
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import curve_fit

class ComputeHydrogenBondFraction:
    """
    compute hydrogen bonding
    """
    def __init__(self, temp, x, y, de_w=None, ds_w=None, de_p=None, ds_p=None, N_list=None, poly_repeats=1, vol_list=None, h_list=[0, 0]):
        """
        fits change in energy and change in entropy for a polymer

        :param temp: ndarray containing simulation temperatures
        :param x: fractions of fraction of polymer hydrogen bonds
        :param y: fractions of fraction of water hydrogen bonds
        :param de_w: change in energy for water-water hydrogen bonds
        :param ds_w: change in entropy for water-water hydrogen bonds
        :param de_p: change in energy for peo-water hydrogen bonds
        :param ds_p: change in entropy for peo-water hydrogen bonds
        :param N_list: list containing number of: polymers, water molecules, cations, anions
        :param poly_repeats: number of polymer repeat units
        :param vol_list: list containing molecular volume of: polymers, water molecules, cations, anions
        :param h_list: list containing nearest neighbors for: cations, anions

        :return: a numpy array
        """

        # initialize parameters
        self.temp = temp
        self.x = x
        self.y = y

        self.kB = 8.31446 # boltzmann's constant (kJ / mol*K)
        self.dE_w, self.dS_w = de_w, ds_w # water change in energy (kT units), water change in entropy (unitless)
        self.de_p, self.ds_p = de_p, ds_p # polymer change in energy (kT units), water change in entropy (unitless)

        Np, Nw, Nc, Na = N_list # number of polymers, number of water molecules, number of salts pairs
        Nr = poly_repeats
        v_p, v_w, v_c, v_a = vol_list # peo repeat unit volume, water molecule volume
        vol_p, vol_w, vol_c, vol_a = v_p*Nr*Np, v_w*Nw, v_c*Nc, v_a*Na # net per molecule type
        V = vol_p + vol_w + vol_c + vol_a # total volume

        self.f_v = v_w/v_p # fraction of water to peo volume

        self.phi_p = vol_p / V # fraction of peo to total volume
        self.phi_w = vol_w / V # fraction of water to total volume

        hc, ha = h_list # cation, anion nearest neihbors
        self.fc, self.fa= hc*(Nc / Nw), ha*(Na / Nw) # cation, anion, fractions

        self.fit_temp = np.append(self.temp, self.temp)
        popt, pcov= curve_fit(self.PolymerSolution, self.fit_temp, np.append(self.x, self.y))
        print("dE_p (kJ/mol), dS_p (unitless):", popt[0], popt[1])

        if self.de_p == None:
            self.de_p, self.ds_p = popt

    def PolymerSolution(self, T, dE_p, dS_p):
        """
        Solve equations determining the fraction of hydrogen bonds given a change in energy and change in entropy

        :param T: temperatures for which to solve the hydrogen bond fractions
        :param dE_p: change in energy for peo-water hydrogen bonds
        :param dS_p: change in entropy for peo-water hydrogen bonds

        :return: a numpy array
        """

        print(dE_p, dS_p, self.fc, self.fa)
        T = np.array(T[:int(len(T)/2)])

        def eqns(val):
            """
            Equations determining the fraction of hydrogen bonds

            :param val: ndarray containing x,p
            :return: equations (ndarray)
            """

            x, y = np.array(val[:len(T)]), np.array(val[len(T):])

            fac = 2 * self.phi_w * (1 - self.fa - y - x * self.f_v * self.phi_p / self.phi_w)

            eqn1 = np.exp(dE_p/(self.kB*T))*np.exp(- dS_p) * (1 - x) * fac
            eqn2 = np.exp(self.dE_w/(self.kB*T))*np.exp(- self.dS_w) * (1 - self.fc - y) * fac

            return np.append(x - eqn1, y - eqn2)


        def solv_eqns(x, y):
            """
            Solution to the equations defining the fraction of hydrogen bonds

            :param x: initial value for fraction of polymer hydrogen bonds
            :param y: initial value for fraction of water hydrogen bonds
            :return: number of hydrogen bonds ( OptimizeResultsObject_ )

            .. _OptimizeResultsObject: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
            """
            sol = fsolve(eqns, np.array([x, y]))

            return sol

        return solv_eqns(self.x, self.y)


    #'''
    def PlotFit(self, size=0):
        """
        plots polymer and water hydrogen fractions with fitted parameters

        :param size: how many points to plot

        :return: None
        """

        x_sim, y_sim = self.x, self.y
        self.x = np.linspace(x_sim[0], x_sim[-1], size)
        self.y = np.linspace(y_sim[0], y_sim[-1], size)
        numer_temp = np.linspace(self.temp[0], self.temp[-1], size)
        numer_temp = np.append(numer_temp, numer_temp)

        mixture = 'PEO+'*len(np.nonzero(self.phi_p!=0)[0]) + 'NaCl+'*len(np.nonzero(self.fc!=0)[0]) + 'H2O'

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', serif='cm10', weight='bold', size=10)
        numer_sol = self.PolymerSolution(numer_temp, self.de_p, self.ds_p)
        plt.plot(numer_temp[:size], numer_sol[:size], label= 'x(T) numerical: '+r'$ \frac{\Delta E_p}{k_B}= $'+ str("{:.2E}".format(self.de_p/self.kB))+' K, '+r'$ \Delta S_p= $'+ str("{:.3}".format(self.ds_p)))
        plt.plot(numer_temp[size:], numer_sol[size:], label=r'$y(T)$'+' numerical: '+r'$\frac{\Delta E_w}{k_B}= $'+ str("{:.2E}".format(self.dE_w/self.kB))+' K, '+r'$ \Delta S_w= $'+ str("{:.3}".format(self.dS_w)))
        plt.plot(self.temp, x_sim, marker='o', ls='None', label='x(T) simulation')
        plt.plot(self.temp, y_sim, marker='o', ls='None', label=r'$y(T)$'+' simulation')
        plt.title(' hydrogen bond fraction, ' + str(mixture))
        plt.ylabel('x(T), y(T) ')
        plt.xlabel('T (K)')
        plt.legend()
        plt.savefig("hb_frac_fit")
        plt.show()


        return None
    #'''

