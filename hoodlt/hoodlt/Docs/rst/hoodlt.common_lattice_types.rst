.. _HOODLTExplainLatticeTypes:

Notes on adding bonds for free energy calculations
==================================================

In order to run a lattice and compute free energies, see :ref:`NcsExample4` and :ref:`SimulationExample5` for example,
it is necessary to add bonds, this is done by a statement such as

.. code-block:: python

    # an empty ConfigurationBuilder object
    builder = ConfigurationBuilder()
    # build the lattice, core objects and ligand objects are passed to ConfigurationBuilder in lists
    builder.build_lattice(list_core_objs, list_lig_objs, lat)

    # add bonds 'CTR-CTR1' for A-B, 'CTR-CTR2' for A-A, 'CTR-CTR3' for B-B
    builder.conf.add_bonds({'CTR-CTR1':(0,1), 'CTR-CTR2':(0,2), 'CTR-CTR3':(1,1)})

The quantity (0,1) means add all degree 1 bonds for particles of type 0 and assign the bond 'CTR-CTR1'.
Degree 1 is nearest neighbors. Similarly (0,2) means add all degree 2 (next to nearest neighbor) and
assign 'CTR-CTR2' and (1,1) means add all nearest neighbor of particles of type 1 and assign 'CTR-CTR3' bonds.

The CTR-CTRn are harmonic bonds defined in the force field. If the default value needs to be changed
it is illustrated, see :ref:`SimulationExample5` by the snippet below

.. code-block:: python

    # here we set the coefficents for both types of CTR-CTR bonds in our system
    # rinit is the equilibrium distance, defaulting to the current separation
    # lamda is a strength scaling factor, defaulting to 1
    sim.initialize('l_00p2', lamda=[0.02, 0.02, 0.02], log_hist='all')

This will change the default value to :math:`\lambda\cdot 0.02` for each of the three bonds. The updated
results will be recorded in the json files. The bond length :math:`r_0` is calculated automatically, as
illustrated below.

Running the following script shows internally how this is done, what lattice points,
what bonds are selected, bond distances, etc..

.. code-block:: python

    import hoodlt.Utils.LatticeNeighbors as Ln
    from hoodlt.Lattices.Lat2.MgZn2_lat import LatMgZn2Base12 as mgzn2

    # create a mgzn2 lattice object with lattice constant 100 and gamma=0.7
    lat = mgzn2(l_value=[2,2,1], a_nn_e=100, gamma=0.7)
    print('created a MgZn2 lattice object with [2,2,1]=4 unit cells')
    print('The lattice constant for this lattice is 100, at gamma = 0.7')

    # instantiate the LatNeighbor class
    lat_neighbor = Ln.LatNeighbor(lat)

    #print the number of lattice points
    print('MgZn2 lattice containing %d particles'%lat_neighbor.l_vec.shape[1])
    print('this is a binary lattice, there are two types (0,1), here is the type of every particle')
    print(lat_neighbor.typ)
    print('Primitive vectors are labelled 0 to 11, so here are the primitive vectors')
    print(lat_neighbor.base)

    print('now we will find all nearest neighbors(degree 1) of point 0, they are', lat_neighbor.neighbor(0,1))
    print('now we will find all next to next to nearest (degree 3) neighbors of point 5, they are', lat_neighbor.neighbor(5,3))
    print('here are the distances for all degres for point 10')
    print(lat_neighbor.u_dist[10])
    print('\n')

    print('this provides all nearest neighbor bonds for particles of type 0')
    print('all nearest neighbors are of type 1')
    l1 = lat_neighbor.bonds_in_lattice(0, 1)
    print(l1)
    print('the distance between these neighbors is', lat_neighbor.u_dist[0][1])
    print('\n')

    print('this provides all nearest neighbor bonds for particles of type 1')
    print('all nearest neighbors are also of type 1')
    l2=lat_neighbor.bonds_in_lattice(1, 1)
    print(l2)
    print('the distance between these neighbors is', lat_neighbor.u_dist[4][1])
    print('\n')

    print('this provides all next to nearest neighbor bonds for particles of type 0')
    print('all nearest neighbors are of type 0')
    l3=lat_neighbor.bonds_in_lattice(0,2)
    print('the distance between these neighbors is', lat_neighbor.u_dist[0][2])
    print(l3)
    print('\n')

    print('from the above all nearest neighbor bonds of particles of type 0 are of type 1')
    print('so (0,1) should be the same as (1,2)')
    l4=lat_neighbor.bonds_in_lattice(1,2)
    print('the distance between these neighbors is', lat_neighbor.u_dist[4][2])
    print(l4)
    print('\n')

