"""
:module: overrp_pot_cutoff
:platform: Unix, Windows
:synopsis: Defines the potential :math:`\\varepsilon(\\frac{\sigma}{r})^p`

    .. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, March2014
"""

from hoodlt.D_matrix.Dmatrix_Mixt_Potential_Simple import DPotential


class OverrpCut(DPotential):

    """Implements the inverse power law potential defined by:

        .. math::
            V(r) = \\varepsilon\\left(\\frac{\sigma}{r}\\right)^p
    """

    def __init__(self, p_e, sig_e, eps_e, cut_e, eps_zero=1e-9):
        """The Constructor

        :param p_e:  p-exponent
        :param sig_e: sigma radius
        :param eps_e: epsilon
        :param cut_e: cut_off
        :param eps_zero: error :math:`|\\vec{R}| < \\varepsilon`
        """

        super(OverrpCut, self).__init__()

        self.p = p_e
        self.sig = sig_e
        self.sig2 = sig_e*sig_e
        self.eps = eps_e
        self.cut = cut_e
        self.epsilon = eps_zero

    def __str__(self):
        """Potential name

        :return: name
        :rtype: str
        """

        nam = 'inverse_power_law ' + 'p=%d' % self.p + 'eps = %1.4f' % self.eps

        return nam

    def der0(self, r):
        """Potential value

        :param r: point
        :return: value
        :rtype: float
        """

        return self.eps/(r/self.sig)**self.p

    def der1(self, r):
        """Potential 1st-derivative value

        :param r: point
        :return: function value
        :rtype: float
        """

        return -(self.eps/self.sig)*self.p/(r/self.sig)**(self.p+1)

    def der2(self, r):
        """Potential 2nd-derivative value

        :param r: point
        :return: function value
        :rtype: numpy.array
        """

        return (self.eps/self.sig2)*self.p*(self.p+1)/(r/self.sig)**(self.p+2)

    def der3(self, r):
        """Potential 2nd-derivative value

        :param r: point
        :return: function value
        :rtype: numpy.array
        """

        return -self.p*(self.p+1)*(self.p+2)*self.eps/(self.sig*self.sig2)/(r/self.sig)**(self.p+3)
