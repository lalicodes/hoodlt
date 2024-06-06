"""
:module: MgZn2_OTM
:platform: Unix, Windows
:synopsis: Defines the clas implementing the OTM :math:`\\mbox{NaZn}_{13}` lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2016
"""

import numpy as np

from scipy.optimize import brenth

from hoodlt.Lattices.Lat2.NaZn13_lat import LatNaZn13Base112
import hoodlt.OTM.OTMobs as ObT


class OTMLatNaZn13Base112(LatNaZn13Base112):
    """ Implementation of the :math:`\\mbox{NaZn}_{13}` lattice
    one hundred twelve particles per unit cell
    """

    def __init__(self, l_value, a_nn_e, gamma, xtol=2.0*1e-12):
        """The constructor

         :param l_value: lattice size
         :param a_nn_e: spacing
         :param gamma: ratio :math:`\\gamma = r_B/r_A , r_B < r_A` of the two particle radii
         :param xtol: tolerance for root finding
         """
        
        self.xtol = xtol
        self.tau = 0.5*(np.sqrt(5)+1)

        self.otm = ObT.OTMObs()
        super(OTMLatNaZn13Base112, self).__init__(l_value, a_nn_e, gamma)

    def otm_observables(self):
        """
        returns the otm observables
        
        :return: The OTM parameters for this lattice
        :rtype: :class:`hoodlt.OTM.OTMobs`
        """

        if self.gam <= self.gamma_crit[2]:
            self.otm.indx = 0
            self.otm.gamma_bar = self.gam
            self.otm.pf = self.pf()

        else:
            def fun(y):
                return self.eqn(y)
                
            self.otm.indx = 1
            self.otm.gamma_bar = brenth(fun, self.gamma_crit[2], self.gam, xtol=self.xtol)
            self.otm.pf = 4.0*np.pi*(1.0+13.0*self.gam**3)/(3.0*(self.a_h(self.otm.gamma_bar))**3)
            self.otm.ratio_b = self.otm.gamma_bar/self.gam

        if self.otm.pf > 1:
            self.otm.indx = -1

        return self.otm

    def a_h(self, y):
        """Helper function, defines the lattice constant in dimensionless numbers
        
        :param y: value of :math:'{\\bar \\gamma}'
        :return: value of the lattice constant
        :rtype: float
        """
        v1 = np.sqrt(self.tau**2+1)
        v2 = np.sqrt(self.tau**2-1)
        return 2.0*(self.gam+y)*(self.tau+v2/(1+self.gam/y))/v1
        
    def eqn(self, y):
        """Helper function, equation defining :math:'{\\bar \\gamma}'
        
        :param y: value of :math:'{\\bar \\gamma}'
        :return: value of the lattice constant
        :rtype: float
        """
        
        c_fac = (1+self.tau)/(4.0*np.sqrt(1+self.tau**2))*self.a_h(y)      
        
        v1 = (self.gam + y)**2/4.0 + 3.0*(self.a_h(y))**2/16.0
        
        return v1-c_fac*(y+self.gam)-0.25*(1+self.gam)**2
