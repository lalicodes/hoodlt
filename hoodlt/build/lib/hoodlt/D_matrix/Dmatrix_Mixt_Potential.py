"""
:module: D_matrix_Mixt_Potential
:platform: Unix, Windows
:synopsis: Defines the abstract class for the lattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, June2014
"""


class DMixtPotential(object):
    """ class defining the potential of a mixture

    """

    def __init__(self):
        self.mrank = None
        self.pot = []
        self.cut = 0.0
        self.epsilon = 0.0

    def __str__(self):
        """Potential name

        :return: potential name
        :rtype: str
        """

        name = ''

        for n in range(self.mrank):
            for m in range(self.mrank):
                name += ' (' + self.pot[n][m].__str__ + ')  '

        return name

    def der0(self, r, sp):
        """The two body potential value

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].der0(r)

    def der1(self, r, sp):
        """The first derivative of the potential

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].der1(r)

    def der2(self, r, sp):
        """The second derivative of the potential

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].der2(r)

    def der3(self, r, sp):
        """The third derivative of the potential

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].der3(r)

    def e1(self, r, sp):
        """The function :math:`e_1(r)`

        :param r: value of r
        :param sp: index of interaction
        :return: value
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].e1(r)

    def e2(self, r, sp):
        """The function :math:`e_2(r)`

        :param r: value of r
        :param sp: index of interaction
        :return: value
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].e2(r)

    def g1(self, r, sp):
        """The function :math:`g_1(r)`

        :param r: value of r
        :param sp: index of interaction
        :return: value
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].g1(r)

    def g2(self, r, sp):
        """The function :math:`g_2(r)`

        :param r: value of r
        :param sp: index of interaction
        :return: value
        :rtype: numpy.array
        """
        return self.pot[sp[0]][sp[1]].g2(r)

    def g2_over_r2(self, r, sp):
        """The function :math:`g_2(r)//r^2`

        :param r: value of r
        :param sp: index of interaction
        :return: value
        :rtype: numpy.array
        """
        return self.pot[sp[0]][sp[1]].g2_over_r2(r)

    def g3(self, r, sp):
        """The function :math:`g_3(r)`

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].g3(r)

    def g3_over_r2(self, r, sp):
        """The function :math:`g_3(r)//r^2`

        :param r: three dimensional numpy array specifying the point
        :param sp: index of interaction
        :return: value of the potential at point r
        :rtype: float
        """

        return self.pot[sp[0]][sp[1]].g3_over_r2(r)
