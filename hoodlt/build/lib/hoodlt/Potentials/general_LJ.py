"""
:module: general_LJ
:platform: Unix, Windows
:synopsis: The Lennard Jones potential :math:`c(p,q)\\varepsilon[(\\sigma/r)^p-(\\sigma/r)^q]`

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
"""

from hoodlt.D_matrix.Dmatrix_Mixt_Potential_Simple import DPotential
from hoodlt.Potentials.overrp_pot_cutoff import OverrpCut


class LJGen(DPotential):
    """Implements the generalized LJ potential defined by

     .. math::
        V_{LJ}(r)=c(p,q)\\varepsilon[(\\sigma/r)^p-(\\sigma/r)^q]
    """

    def __init__(self, p_e, q_e, sig_e, eps_e, cut_e, eps_zero=1e-10):
        """The construtor

        :param p_e: power of the repulsive part
        :param q_e: power of the attractive part
        :param sig_e: sigma radius
        :param eps_e: epsilon coefficient
        :param cut_e: cut_off
        :param eps_zero: defines :math:`|\\vec{R}| < \\varepsilon`
        """

        super(LJGen, self).__init__()

        self.cut = cut_e

        self.p = p_e
        self.q = q_e
        self.sigma = sig_e
        self.e = eps_e

        p_f = float(p_e)
        q_f = float(q_e)
        self.const = eps_e*p_f/((p_f-q_f)*(q_f/p_f)**(q_f/(p_f-q_f)))

        self.pot_p = OverrpCut(p_e, sig_e, 1.0, cut_e, eps_zero)
        self.pot_q = OverrpCut(q_e, sig_e, 1.0, cut_e, eps_zero)

    def __str__(self):
        """Potential name

        :return: name
        :rtype: str
        """

        nam = 'Generalized LJ ' + 'p=%f' % self.p + 'q = %f' % self.q + 'sigma = %f' % self.sigma + 'eps = %f' % self.e

        return nam

    def der0(self, r):
        """Potential value

        :param r: point
        :return: value
        :rtype: float
        """

        return self.const*(self.pot_p.der0(r)-self.pot_q.der0(r))

    def der1(self, r):
        """Potential 1st-derivative value

        :param r: point
        :return: value
        :rtype: float
        """
        return self.const*(self.pot_p.der1(r)-self.pot_q.der1(r))

    def der2(self, r):
        """Potential 2nd-derivative value

        :param r: point
        :return: value
        :rtype: numpy.array
        """
        return self.const*(self.pot_p.der2(r)-self.pot_q.der2(r))

    def der3(self, r):
        """Potential value pressure

        :param r: point
        :return: value
        :rtype: float
        """
        return self.const*(self.pot_p.der3(r)-self.pot_q.der3(r))
