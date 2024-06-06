"""
:module: overrp_pot_cutoff_vanish_Mixt
:platform: Unix, Windows
:synopsis: Defines the potential of a mixture of inverse power laws and zero potentials

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, July2014
"""

from hoodlt.D_matrix.Dmatrix_Mixt_Potential import DMixtPotential
from hoodlt.Potentials.overrp_pot_cutoff import OverrpCut


class OverrpCutVanishMixt(DMixtPotential):
    """Defines particles that interact with inverse power laws and zero potentials. This class is suitable for unit
    tests
    
    """

    def __init__(self, p, ind, sig, eps, cut_off, eps_zero=1e-10):
        """The constructor

        :param p: exponent of the potential
        :param ind: index of the pair potential that is not zero
        :param sig: is a list whose n elements are :math:`\\sigma_1,\\cdots,\\sigma_a,\\cdots,\sigma_n`
        :param eps: s a list whose :math:`n^2` elements are the :math:`\\varepsilon_{ab}` coefficients.
        :param cut_off: cut-off value
        :param eps_zero: defines :math:`|\\vec{R}| < \\varepsilon`
        :return:
        """

        super(OverrpCutVanishMixt, self).__init__()

        self.epsilon = eps_zero
        self.mrank = len(sig)
        self.cut = cut_off
        self.epsilon = eps_zero

        for n in range(self.mrank):
            self.pot.append([])
            for m in range(self.mrank):
                if (n == ind[0]) and (m == ind[1]):
                    vp = OverrpCut(p, 0.5*(sig[n]+sig[m]), eps[n, m], cut_off, eps_zero)
                else:
                    vp = OverrpCut(p, 0.5*(sig[n]+sig[m]), 0.0, cut_off, eps_zero)
                self.pot[n].append(vp)
