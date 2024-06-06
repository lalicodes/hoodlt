"""
:module: SphereInfo
:platform: Unix, Windows
:synopsis: Defines the class for Spherical particle Information

.. moduleauthor:: Nathan Horst <nhorst@iastate.edu>, June 2017
"""

from __future__ import division


class SphereInfo(object):
    """Initialization information for a general spherical nanoparticle
    """

    def __init__(self, core_radius, graft_radius, core_num, graft_num):

        self.core_radius = core_radius
        self.graft_radius = graft_radius
        self.core_num = core_num
        self.graft_num = graft_num
