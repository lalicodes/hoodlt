.. _LatticesUtilities:

Analysis of lattices
====================

Lattice utilities allow to find Voronoi cells, compute displacements and provide different
ways to query lattices about certain lattice sites and obtain, for example, the number of
nearest neighbors, next to nearest neighbors, etc..

A point in a lattice maybe characterized by five integers: :math:`(i_b, j_b, l_1, l_2, l_3)`

.. math::

      {\bm R}_{(i_b,j_b,l_1,l_2,l_3)}={\bm v}_{(i_b,j_b)}+l_1 {\bm a}_1 + l_2 {\bm a}_2 + l_3 {\bm a}_3

where :math:`(i_b,j_b)` parameterizes the bases within the unit cell
and :math:`(l_1,l_2,l_3)` the translations of the unit cell. For example, in a binary lattice
:math:`i_b=0,1` and :math:`j_b = 0 \cdots n_{i_b}-1` where :math:`n_{i_b}` is the number of
basis of type :math:`i_b`. There are :math:`n_{total}=n_{basis} L_1 L_2 L_3` lattice sites,
where :math:`n_{basis}=\sum_{i_b=0} n_{i_b}` is the number of basis within the unit cell.

    Alternatively, each point may be characterized by a singe integer :math:`i_{val} , 0 \le i_{val} \le n_{total}-1`

.. math::

    i_{val}(i_b,j_b,l_1,l_2,l_3) = j_b+\sum_{i=0}^{i=i_b-1} n_i + n_{basis}(l_1 + L_1 l_2 + L_1 L_2 l_3)

with :math:`l_i = 0 \cdots L_i-1`. Alternatively, given :math:`i_{val}` it is possible to
obtain :math:`(i_b,j_b,l_1,l_2,l_3)` recursively from the formula:

.. math::
    \begin{array}{c c c}
    l_3 &=& \frac{i_{val}}{n_{basis} L_1 L_2} \\
    l_2 &=& \frac{i_{val} - n_{basis} L_1 L_2 l_3 }{n_{basis} L_1} \\
    l_1 &=& \frac{i_{val} - n_{basis} (L_1 l_2  + L_1 L_2 l_3)}{n_{basis}} \\
    j_b+\sum_{i=0}^{i=i_b-1} n_i  &=& i_{val}-n_{basis}(l_1+ L_1 l_2 + L_1 L_2 l_3)
    \end{array}

**Example:**

Let us identify the neighbors of two sites in :math:`\mbox{Mg}\mbox{Zn}_2` structure, output
what the distance is and the coordinates of all points in the different formats

.. code-block:: python

    import numpy as np

    import hoodlt.Lattices.Lat2.MgZn2_lat as MgZn2_lat
    import hoodlt.Utils.LatticeNeighbors as LatN
    import hoodlt.Utils.LatticePoints as LatP

    # define gamma
    gam = 0.75

    # create the MgZn2 BNSL
    lat = MgZn2_lat.LatMgZn2Base12(3, 1.0, gam)

    # LatNeighbor object
    latn = LatN.LatNeighbor(lat)

    # LatPoints object
    latp = LatP.LatPoints(lat)

    # basis
    ind_a = 3
    ind_b = 7
    # choose this unit cell
    l_1 = 0
    l_2 = 0
    l_3 = 0

    # type of nearest neighbors 1=nn, 2=nnn, etc..
    ind_neigh = 2

    print('A')
    # get all neighbors of point (0, ind_a, l_1, l_2, l_3)
    ind_p = latp.i_val(np.array([0, ind_a, l_1, l_2, l_3]))
    # get all neighbors of index ind_neigh
    neigh_five = latn.neighbor(ind_p, ind_neigh)
    # print the result
    print(ind_a, neigh_five)
    # print the distance between the actual point and its neighbors
    print('distance %1.4f'%latn.u_dist[ind_p][ind_neigh])
    # print the 5 coordinate of each one of the neighbors
    for ind_p in neigh_five:
        crd = latp.coordinate_five(ind_p)
        print(ind_p, crd)

    print('\n')

    print('B')
    # get all neighbors of point (0, ind_b, l_1, l_2, l_3)
    ind_p = latp.i_val(np.array([1, ind_b, l_1, l_2, l_3]))
    # get all nearest neighbors
    neigh_five = latn.neighbor(ind_p, ind_neigh)
    # print the result
    print(ind_b, neigh_five)
    # print the distance between the actual point and its neighbors
    print('distance %1.4f'%latn.u_dist[ind_p][ind_neigh])
    # print the 5 coordinate of each one of the neighbors
    for ind_p in neigh_five:
        print(ind_p, latp.coordinate_five(ind_p))



LatticePoints
-------------

.. automodule:: hoodlt.Utils.LatticePoints
    :members:
    :undoc-members:
    :show-inheritance:

LatticeCoords
-------------

.. automodule:: hoodlt.Utils.LatticeCoords
    :members:
    :undoc-members:
    :show-inheritance:

Lattice_Geometry
----------------

.. automodule:: hoodlt.Utils.Lattice_Geometry
    :members:
    :undoc-members:
    :show-inheritance:

LatticeNeighbors
----------------

.. automodule:: hoodlt.Utils.LatticeNeighbors
    :members:
    :undoc-members:
    :show-inheritance:


Lattice_UnitCell
----------------

.. automodule:: hoodlt.Utils.Lattice_UnitCell
    :members:
    :undoc-members:
    :show-inheritance:

GeomClasses
-----------

.. automodule:: hoodlt.Utils.GeomClasses
    :members:
    :undoc-members:
    :show-inheritance:

voro_hoodlt
-----------

.. automodule:: hoodlt.Utils.voro_hoodlt
    :members:
    :undoc-members:
    :show-inheritance:

qhull_hoodlt
------------

.. automodule:: hoodlt.Utils.qhull_hoodlt
       :members:
    :undoc-members:
    :show-inheritance:

Displacements
-------------

.. automodule:: hoodlt.Utils.Displacements
    :members:
    :undoc-members:
    :show-inheritance:

