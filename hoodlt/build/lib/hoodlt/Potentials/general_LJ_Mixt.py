"""
:module: general_LJ_Mixt
:platform: Unix, Windows
:synopsis: Defines the potential of a mixture of Lennard Jones potentials

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November2018
"""

from hoodlt.D_matrix.Dmatrix_Mixt_Potential import DMixtPotential
from hoodlt.Potentials.general_LJ import LJGen


class LJMixt(DMixtPotential):
    """ Defines class OverrpCutMixt, defined by

    .. math::
        V_{a,b}(r)=\\varepsilon_{ab}\\left(\\frac{\sigma_{ab}}{r}\\right)^p

    .. math::
       \\sigma_{ab}=(\\sigma_{a}+\\sigma_{b})/2

    :math:`\\varepsilon_{ab}` is a nxn matrix where n is the number of species
    """

    def __init__(self, p, q, sig, eps, cut_off, eps_zero=1e-10):
        """The constructor

        :param p: exponent of the repulsive part of the potential
        :param q: exponent of the attractive part of the potnetial
        :param sig: is a list whose n elements are :math:`\\sigma_1,\\cdots,\\sigma_a,\\cdots,\sigma_n`
        :param eps: is a list whose :math:`n^2` elements are the :math:`\\varepsilon_{ab}` coefficients.
        :param cut_off: cut-off value
        :param eps_zero: defines :math:`|\\vec{R}| < \\varepsilon`
        :return:
        """

        super(LJMixt, self).__init__()

        self.epsilon = eps_zero

        self.mrank = len(sig)
        self.cut = cut_off

        for n in range(self.mrank):
            self.pot.append([])
            for m in range(self.mrank):
                self.pot[n].append(LJGen(p, q, 0.5*(sig[n]+sig[m]), eps[n, m], cut_off, eps_zero))
