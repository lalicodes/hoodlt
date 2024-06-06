"""
:module: pf_sphere_cube
:platform: Unix, Windows
:synopsis: Function computing the packing fraction for a A-sphere/ B-cube system

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np


def pf_s_c(self):
    """
        Corrected packing fraction
    """
    return self.pf()*(self.typ[0]+6*self.typ[1]*self.gam**3/np.pi)/(self.typ[0]+self.typ[1]*self.gam**3)
