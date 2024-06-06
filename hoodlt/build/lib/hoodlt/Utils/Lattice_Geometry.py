"""
:module: Lattice_Geometry
:platform: Unix, Windows
:synopsis: Functions to calculate different geometric properties of a lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, April2015
"""
import pickle


def naming_convention(v_data, ind):
    """name of an output file for point ind in a VoroData object

    :param v_data: a VoroData object
    :param ind: index of base point
    :return: file name
    :rtype: str
    """

    lat = v_data.u_c.lat
    f_name = lat.name()[0]
    f_name += '_' + str(v_data.u_c.u_cell[ind, 1]) + ',' + str(v_data.u_c.u_cell[ind, 2]) + '_'

    # name file according to diameter not radius
    l_name = [2.0*rad for rad in lat.radius]

    f_name += '_'.join('%1.4f' % v for v in l_name)

    return f_name


def draw_facets(v_data, ind):
        """Returns a structure that can be used to draw the unit cell

        :param v_data: a VoroData object
        :param ind: index of the base point
        """

        fid = open(naming_convention(v_data, ind), 'w')

        pickle.dump(v_data.ws_simple(ind), fid)
