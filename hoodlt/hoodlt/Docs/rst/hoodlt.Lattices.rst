.. _ListLattices:

Lattices
********

The Lattice class is fundamental to HOODLT. It contains the necessary information to build
all points in a lattice, but does not store any lattice positions.

The Lattice class is designed to be lightweight and contain the minimal information.
Utilities that provide more detailed and selective information are given in:
Utils, see Analysis of Lattices in :ref:`ListUtilities`

Example
-------

Let us create a :math:`\mbox{cubAB}_{13}` lattice object and show all the functionalities
that the class provides. This is done with the snippet.

.. code-block:: python

        import hoodlt.Lattices.Lat2.cub_AB13_lat as CuBaB

        l_box = 3
        # define a box with 3 x 3 x 3 unit cells

        a_nn = 1 
        # minimum separation between A-particles 

        gamma = 0.8
        # ratio of the diameter of B and A particle

        cub = CuBaB.LatCubAB13Base14(l_box, a_nn, gamma)
        # create lattice object

        print(cub.a_val(), cub.pf(), cub.num_pnts())
        # prints a_nn, packing fraction, total number of particles

        print(cub.l_box(), cub.vol_unit_cell())
        # prints box size (in the same units as a_nn), volume unit cell

        print(cub.l, cub.typ, cub.radius)
        # prints the number of unit size, particle types and the radius        

        print(cub.get_a(1))
        # prints the 2nd primitive vector as a 3d numpy array 

        print(cub.get_v([1, 0]))
        # prints the first vector base corresponding to a B-particle (type 1)
        print(cub.get_v([0, 1]))
        # prints the second vector base corresponding to a A-particle (type 0)

        cub_name = cub.name()
        print(cub_name)
        # prints a list that contains, in this order:
        # Lattice name, LaTeX name, space group, 
        # Strukturberich and its assigned color.


Subpackages
-----------

.. toctree::

    hoodlt.Lattices.Lat1
    hoodlt.Lattices.Lat2
