"""
:module: PlotLatticePotential
:platform: Unix, Windows
:synopsis: Defines the class used to plot PMF results

.. moduleauthor: Xun Zha <xzha@iastate.edu> July 2019
.. history:
..                Xun Zha <xzha@iastate.edu> July 2021
..                  - bug fix; added new method animation_free_energy()
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.animation import ArtistAnimation


class PlotLatticePotential:
    """
    This class plots lattice potentials and other thermodynamic functions
    """
    def __init__(self, clp_class):
        """Initializer

        :param clp_class: ComputeFreeEnergyLattice class object
        """

        # initialize
        self.clp = clp_class

        self.lat_typ = self.clp.lat_typ
        self.a_val = self.clp.a_val
        self.a_nn = self.clp.get_a_nn()

        # use tex in text
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif', serif='cm10', weight='bold', size=14)
        plt.rc('legend', fontsize=12)

    def plot_pressure(self):
        """Plots the pressure as a function of a_nn

        :return:
        """

        # number of lattice points
        num_pnts = sum(self.clp.f_object0.lat.num_pnts())

        # pressure
        pressure = self.clp.get_pressure()

        # average volume occupied by one NC
        volume_per_nc0 = self.clp.f_object0.lat.vol_unit_cell() / np.sum(self.clp.f_object0.lat.typ)
        volume_per_nc = np.zeros(len(self.a_val))
        for i in range(len(volume_per_nc)):
            volume_per_nc[i] = volume_per_nc0 * np.power(self.clp.a_ratio[i], 3)
        fit_volume_per_nc = np.linspace(volume_per_nc[0], volume_per_nc[-1], 2000)

        # the spring contribution to pressure
        if self.lat_typ == 'snsl':
            spring_virial = self.clp.virial_data
        else:
            spring_virial = self.clp.virial_data[0]
        adjust_pressure = np.mean(spring_virial, axis=1) / 3 / volume_per_nc / num_pnts

        plt.figure()
        plt.title('Pressure')
        contact_name = self._get_contact_name()
        plt.xlabel(contact_name + r'$a_{nn}$ ($\mbox{\AA}$)')
        plt.ylabel('Pressure (MPa)')
        plt.plot(self.a_nn, (pressure + adjust_pressure) * self.clp.pressure2mpa, 'o', label='Pressure all')
        plt.plot(self.a_nn, pressure * self.clp.pressure2mpa, '*', label='Pressure subtracting harmonic')
        if self.clp.pp is not None:
            fit_a = np.linspace(self.a_nn[0], self.a_nn[-1], len(self.clp.pp))
            pair_pressure = -np.gradient(self.clp.pp, fit_volume_per_nc) * self.clp.pressure2mpa
            plt.plot(fit_a, pair_pressure, c='k', alpha=.4, label='pair result')
        plt.legend(loc='lower right')
        plt.tight_layout()

    def plot_virial(self):
        """Plots the virial and energy contribution from the harmonic bonds

        :return:
        """

        # spring contribution
        spring_virial = np.mean(self.clp.virial_data, axis=-1)
        spring_energy = np.mean(self.clp.energy_data, axis=-1)

        # plot virial and spring energy
        plt.figure()
        plt.title('Harmonic bond energy and virials per NC')
        contact_name = self._get_contact_name()
        plt.xlabel(contact_name + r'a$_{nn}$ (\AA)')
        plt.ylabel('k$_B$T')
        if self.lat_typ == 'snsl':
            plt.plot(self.a_nn, spring_virial, 'o', label='spring virial')
            plt.plot(self.a_nn, spring_energy, 'ro', label='spring energy')
        else:
            plt.plot(self.a_nn, spring_virial[0], 'o', label='spring virial')
            plt.plot(self.a_nn, spring_energy[0], 'ro', label='spring energy')
            plt.plot(self.a_nn, spring_virial[1], 'o', label='spring virial (type A)')
            plt.plot(self.a_nn, spring_energy[1], 'ro', label='spring energy (type A)')
            plt.plot(self.a_nn, spring_virial[2], 'o', label='spring virial (type B)')
            plt.plot(self.a_nn, spring_energy[2], 'ro', label='spring energy (type B)')
        plt.legend()
        plt.tight_layout()

    def plot_free_energy(self, file_n, title_n, opm=None, otm=None):
        """plots the pressure and the free energy as a function of a_nn

        :param file_n: filename of the output plots
        :param title_n: title of the plots
        :param opm: OPM value
        :param otm: OTM value for binary systems
        :return:
        """

        # plot the free energy
        fig = self._plot_free_energy(title_n, opm, otm)

        # save the plot in pdf format
        pp = PdfPages(file_n+'_free_energy.pdf')
        pp.savefig(fig)
        pp.close()

    def animation_free_energy(self, indices, file_n, title_n, opm=None, otm=None, fps=10, dpi=200,
                              text_x=None, text_y=None):
        """plots the pressure and the free energy as a function of a_nn, create an animation with
        the pressure data point highlighted at a given a_nn of a given frame

        :param indices: indices of the highlighted a_nn's
        :param file_n: filename of the output gif
        :param title_n: title of the plots
        :param opm: OPM value
        :param otm: OTM value for binary systems
        :param fps: frames per second
        :param dpi: dots per linear inch
        :param text_x: position in (x, y) for the highlighted text
        :param text_y: position in (x, y) for the highlighted text
        :return:
        """

        # plot the free energy
        fig = self._plot_free_energy(title_n, opm, otm)
        (ax1, ax2) = fig.axes

        # data points in the animation
        x = self.a_nn[indices]
        y = self.clp.get_pressure()[indices] * self.clp.pressure2mpa  # convert pressure unit to MPa

        # color of the highlighted data points
        color = 'red'

        # text for the highlighted data point
        if text_x is None:
            text_x = np.amax(x)
        if text_y is None:
            text_y = np.amin(y)/3.
        text_prefix = r'a$_{nn}$='
        text_suffix = r'$\mbox{\AA}$'

        # create the artist list
        imgs = []
        for ind in range(len(x)):
            im1, = ax1.plot(x[ind], y[ind], 'o', c=color)
            im2 = ax1.axvline(x=x[ind], alpha=.5, c=color)
            im3 = ax2.axvline(x=x[ind], alpha=.5, c=color)
            im4 = ax1.text(text_x, text_y, text_prefix+str(np.around(x[ind], 1))+text_suffix, fontsize=20, c=color)
            imgs.append([im1, im2, im3, im4])

        # save the animation
        ani = ArtistAnimation(fig, imgs, blit=True)
        ani.save(file_n+'.gif', writer='pillow', fps=fps, dpi=dpi)

    def _plot_free_energy(self, title_n, opm=None, otm=None):
        """plots the pressure and the free energy as a function of a_nn

        :param title_n: title of the plots
        :param opm: OPM value
        :param otm: OTM value for binary systems
        :return:
        """

        # pressure
        pressure = self.clp.get_pressure()
        # free energy
        fit_free_energy = self.clp.get_fit_free_energy()

        # fit a_nn
        fit_a = np.linspace(self.a_nn[0], self.a_nn[-1], len(fit_free_energy))

        # print information about minimum point of free energy
        min_arg = np.argmin(fit_free_energy)
        print('Minimum free energy %6.2f kBT, at a_nn %5.2f A' % (fit_free_energy[min_arg], fit_a[min_arg]))
        # print information about pair potential approximation of free energy
        if self.clp.pp is not None:
            print('Pair potential sum result is %f kBT' % (self.clp.pp[min_arg]))

        # create canvas
        contact_name = self._get_contact_name()
        f0, (ax1, ax2) = plt.subplots(nrows=2)
        ax1.set_title(title_n)
        ax1.set_ylabel('Pressure (MPa)')
        ax2.set_xlabel(contact_name + r'a$_{nn}$ ($\mbox{\AA)}$')
        ax2.set_ylabel('Energy (k$_B$T)')

        # plot pressure
        pressure = pressure * self.clp.pressure2mpa  # convert pressure unit to MPa
        ax1.plot(self.a_nn, pressure, 'o', label='Pressure', c='C0')
        fit_pressure = interpolate.pchip_interpolate(self.a_nn, pressure, fit_a)
        ax1.plot(fit_a, fit_pressure, c='C0')
        ax1.axhline(y=0, color='k', linestyle='--')
        ax1.axvline(x=fit_a[min_arg], c='C0', ls='--', label='Free energy minimum')
        ax1.set_ylim(min(pressure)*1.1, 5)

        # plot Helmholtz free energy
        ax2.plot(fit_a, fit_free_energy, label='Free energy')
        ax2.axhline(y=0, color='k', linestyle='--')
        ax2.axvline(x=fit_a[min_arg], c='C0', ls='--', label='Free energy minimum')

        # pair potential approximation
        if self.clp.pp is not None:
            ax2.plot(fit_a, self.clp.pp, label='Sum of pair potential')
            ax2.plot(fit_a, (fit_free_energy - self.clp.pp), label='Many body effects')

        # predicted values of models for a_min
        left = self.a_nn[0]
        if opm is not None:  # OPM prediction
            ax2.axvline(x=opm, c='k', alpha=0.5, ls='--', label='OPM')
            left = np.amin([opm, left])
        if otm is not None:  # OTM prediction
            ax2.axvline(x=otm, c='r', ls='--', label='OTM')
            left = np.amin([otm, left])
        # adjust the plot
        right = fit_a[np.argmin(np.abs(interpolate.pchip_interpolate(self.a_nn, pressure, fit_a)))]
        ax1.set_xlim(left - 1, right + 1)
        ax2.set_xlim(left - 1, right + 1)
        ax1.legend(loc='lower right', frameon=False)
        ax2.legend(loc='lower right', frameon=False)
        plt.tight_layout()

        return f0

    def _get_contact_name(self):
        """Contact name shown on the plot

        :param contact: String, "AA", "BB" or "AB" for binary systems
        :return: contact name
        """

        if self.lat_typ == 'snsl':
            return ''

        dict_name = {'AA': 'A-NC ', 'AB': 'A-B ', 'BB': 'B-NC '}
        return dict_name[self.clp.contact]
