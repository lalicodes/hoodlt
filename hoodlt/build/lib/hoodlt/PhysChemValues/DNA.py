"""
:module: DNA
:platform: Unix, Windows
:synopsis: returns general properties for DNA in superlattices

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2020
"""

import numpy as np
import numpy.linalg as la


class DNA(object):
    """
    Defines the class defining the parameters for DNA: those are taken from
    Nanoparticle Superlattice Engineering with DNA
    R. Macfarlane et al.
    Science 334 (2011) 204-208.
    """

    def __init__(self, num_basis, r_core):
        """
        Defines the constructor

        :param num_basis: number of base pairs
        :param r_core: radius of the nanoparticle core

        """

        # num_basis includes A10 spacer, recognizing sequence, unpaired flexor basis and sticky end
        self.num_basis = num_basis
        self.r_core = r_core

        # hydrodynamic radius of a DNA nanoparticle: the 0.4 accounts for the propyl-thiol moiety
        self.hydro_radius_thiol = r_core + 0.34*num_basis + 0.4


def superlattice_ab(dna1, dna2, l_sticky_end):
    """
    Computes the lattice constant for two dna objects: those are taken from
    Nanoparticle Superlattice Engineering with DNA
    R. Macfarlane et al.
    Science 334 (2011) 204-208.

    :param dna1: DNA object
    :param dna2: DNA object
    :param l_sticky_end: length of the complementary sticky ends
    :return : a list with the lattice constant, gamma, overlap parameter, fractions for A and B and the DNA objects
    :rtype: list
    """

    lattice_constant = dna1.r_core + dna2.r_core + 0.255*(dna1.num_basis + dna2.num_basis - l_sticky_end) + 0.8

    if dna1.hydro_radius_thiol > dna2.hydro_radius_thiol:
        dna_a = dna1
        dna_b = dna2
    else:
        dna_a = dna2
        dna_b = dna1

    # gamma value
    gam = dna_b.hydro_radius_thiol/dna_a.hydro_radius_thiol

    # compute overlap parameter = :math:`\\frac{lattice_constant}{r_A+r_B}`
    ovr = lattice_constant/(dna_a.hydro_radius_thiol+dna_b.hydro_radius_thiol)

    # return the fraction overlap for particles A and B
    f_a =  1-(1-gam+ovr**2*(1+gam))/(2*ovr)
    f_b =  1-(-1+gam +ovr**2*(1+gam))/(2*gam*ovr)

    # return these values as a list
    return [lattice_constant, gam, ovr, f_a, f_b, dna_a, dna_b]
