.. _NcsExample4:

Example 4: Lattice and Lattice from Reinitialized Particles
===========================================================

All of the configurations built in :ref:`NcsExample1` through :ref:`NcsExample3` dealt with configurations of
nanoparticles which did not use periodic boundary conditions; lattices are not this simple. To build a lattice, we
will be using the function :func:`hoodlt.Data.Modelconfigurations.ConfigurationBuilder.build_lattice()`,
while utilizing a lattice object :mod:`hoodlt.Data.Modelconfigurations.LatticeFunctionalizedConfiguration`.

Often times when building lattices at a relatively small lattice constant, it is necessary to compress the nanoparticles
that you use to build the lattice so the positions of the ligand chains do not overlap and the simulation can start
successfully, and then build the lattice with reinitialized nanoparticles using the function
:func:`hooldt.Data.Modelconfigurations.ConfigurationBuilder.build_lattice_from_reinit()`, as shown below in examples.


1. Build single-component lattices.

.. code-block:: python

    import numpy as np
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand
    from hoodlt.Lattices.Lat1.fcc_lat_Mixt import LatFccBase4 as fcc
    from hoodlt.Data.ProcessConfigurations.Squeeze import Squeeze

    # name of the forcefield
    forcefield = 'dry-ncs'

    # build a fcc lattice from core object and ligand object
    # fcc lattice object
    lat = fcc(l_value=2, a_nn_e=100)

    # core object
    core = TO201(forcefield)

    # ligand object
    lig = HydrocarbonLigand(11, forcefield)

    # an Empty ConfigurationBuilder object
    builder = ConfigurationBuilder()

    # build the lattice, core object and ligand objects are passed to ConfigurationBuilder in lists
    builder.build_lattice([core], [[lig]*core.graft_num], lat)

    # add bonds 'CTR-CTR1' for nearest neighbors and 'CTR-CTR2' for next nearest neighbors
    builder.conf.add_bonds({'CTR-CTR1':(0,1), 'CTR-CTR2':(0,2)})

    # save configuration to a gsd file and a pickle file
    save_config(builder.conf)

The main issue here is regarding the bonds in :ref:`HOODLTExplainLatticeTypes` a detailed explanation on what
and how bonds  (through add_bonds) are added is provided.

**Views:**

.. raw:: html

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
    <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>

    <style>
      .carousel-caption {
          background-color: rgba(0, 0, 0, 0.5);
          padding: 10px;
      }
      .carousel {
          max-width: 500px;
          margin: auto;
      }
      .carousel-inner img {
          max-width: 100%;
          height: auto;
      }
    </style>

    <div id="carouselExampleIndicators" class="carousel slide" data-ride="carousel">
      <ol class="carousel-indicators">
        <li data-target="#carouselExampleIndicators" data-slide-to="0" class="active"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="1"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="2"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="3"></li>
      </ol>
      <div class="carousel-inner">
        <div class="carousel-item active">
          <img class="d-block w-100" src="_static/latticefront.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeleft.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeperspective.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticetop.png" alt="Top view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Top View</h5>
          </div>
        </div>
      </div>
      <a class="carousel-control-prev" href="#carouselExampleIndicators" role="button" data-slide="prev">
        <span class="carousel-control-prev-icon" aria-hidden="true"></span>
        <span class="sr-only">Previous</span>
      </a>
      <a class="carousel-control-next" href="#carouselExampleIndicators" role="button" data-slide="next">
        <span class="carousel-control-next-icon" aria-hidden="true"></span>
        <span class="sr-only">Next</span>
      </a>
    </div>

2. Build single-component lattices by squeezing ligands first.

We often need to place the nanoparticles close. For that purpose, we use
the function :py:class:`hoodlt.Data.ProcessConfigurations.Squeeze`
whose role is to squeeze the nanoparticles to a given radius so that they can fit into the
lattice. Here is how it works:

.. code-block:: python

    import numpy as np
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand
    from hoodlt.Lattices.Lat1.fcc_lat_Mixt import LatFccBase4 as fcc
    from hoodlt.Data.ProcessConfigurations.Squeeze import Squeeze

    # name of the forcefield
    forcefield = 'dry-ncs'

    # case 2
    # build a fcc lattice from reinitialized FunctionalizedParticle object, (a_nn can be smaller)
    # fcc lattice object
    lat = fcc(l_value=2, a_nn_e=50)

    # build a single NC (a FunctionalizedParticle object)
    # an empty ConfigurationBuilder object for the single NC
    builder = ConfigurationBuilder()
    # core object
    core = TO201(forcefield)
    # ligand object
    lig = HydrocarbonLigand(11, forcefield)
    # build a single NC at position [0, 0, 0]
    builder.add_nc(core, [lig]*core.graft_num, [0, 0, 0])
    # get the FunctionalizedParticle object
    nc = builder.conf.particles[0]

    # squeeze the single NC, FunctionalizedParticle object and radius, shape are passed to Squeeze in lists, nsteps is
    # number of steps to squeeze NC into given radius, a list of FunctionalizedParticle object(s) is returned
    list_nc = Squeeze(forcefield, list_nc=[nc], list_radius=[16]).squeeze()

    # build fcc from reinitialized FunctionalizedParticle object
    # an empty ConfigurationBuilder object for the fcc lattice
    builder = ConfigurationBuilder()
    # FunctionalizedParticle object are passed to builder in list
    builder.build_lattice_from_reinit(list_nc, lat)
    # add harmonic bonds between NCs
    builder.conf.add_bonds({'CTR-CTR1':(0,1), 'CTR-CTR2':(0,2)})

    # save configuration to a gsd file and a pickle file
    save_config(builder.conf)

You should pay attention to the files:

.. code-block:: bash

    cAu201S-Hydrocarbon-n11_cFcc_l2_a50_ffDry-Ncs_uAngAmuEv_bonds.json
    cAu201S-Hydrocarbon-n11_cFcc_l2_a100_ffDry-Ncs_uAngAmuEv_bonds.json

they contain all the information about the bonds CTR-CTR1 and CTR-CTR2 created. Those will be
used both for running simulations as well as to analyze the data.

**Views:**

.. raw:: html

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
    <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>

    <style>
      .carousel-caption {
          background-color: rgba(0, 0, 0, 0.5);
          padding: 10px;
      }
      .carousel {
          max-width: 500px;
          margin: auto;
      }
      .carousel-inner img {
          max-width: 100%;
          height: auto;
      }
    </style>

    <div id="carouselExampleIndicators" class="carousel slide" data-ride="carousel">
      <ol class="carousel-indicators">
        <li data-target="#carouselExampleIndicators" data-slide-to="0" class="active"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="1"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="2"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="3"></li>
      </ol>
      <div class="carousel-inner">
        <div class="carousel-item active">
          <img class="d-block w-100" src="_static/latticefront2.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeleft2.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeperspective2.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticetop2.png" alt="Top view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Top View</h5>
          </div>
        </div>
      </div>
      <a class="carousel-control-prev" href="#carouselExampleIndicators" role="button" data-slide="prev">
        <span class="carousel-control-prev-icon" aria-hidden="true"></span>
        <span class="sr-only">Previous</span>
      </a>
      <a class="carousel-control-next" href="#carouselExampleIndicators" role="button" data-slide="next">
        <span class="carousel-control-next-icon" aria-hidden="true"></span>
        <span class="sr-only">Next</span>
      </a>
    </div>


3. Build binary lattices.

.. code-block:: python

    import numpy as np
    import numpy.linalg as la
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand
    from hoodlt.Lattices.Lat2.MgZn2_lat import LatMgZn2Base12 as mgzn2
    from hoodlt.Data.ProcessConfigurations.Squeeze import Squeeze

    # name of the forcefield
    forcefield = 'dry-ncs'

    # build a binary MgZn2 lattice from core objects and ligand objects
    # MgZn2 lattice object
    lat = mgzn2(l_value=2, a_nn_e=100, gamma=36.77/53.29)

    # core objects
    core_a = TO201(forcefield)
    core_b = TO201(forcefield)
    list_core_objs = [core_a, core_b]

    # ligand objects
    lig_a = HydrocarbonLigand(11, forcefield)
    lig_b = HydrocarbonLigand(11, forcefield)
    list_lig_objs = [[lig_a]*core_a.graft_num, [lig_b]*core_b.graft_num]
    # an empty ConfigurationBuilder object
    builder = ConfigurationBuilder()
    # build the lattice, core objects and ligand objects are passed to ConfigurationBuilder in lists
    builder.build_lattice(list_core_objs, list_lig_objs, lat)

    # add bonds 'CTR-CTR1' for A-B, 'CTR-CTR2' for A-A, 'CTR-CTR3' for B-B
    builder.conf.add_bonds({'CTR-CTR1':(0,1), 'CTR-CTR2':(0,2), 'CTR-CTR3':(1,1)})

    # save configuration to a gsd file and a pickle file
    save_config(builder.conf)


A list of all available lattices and lattice objects is at :ref:`ListLattices`.

Note that here the json files contain three bonds.

**Views**

.. raw:: html

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
    <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>

    <style>
      .carousel-caption {
          background-color: rgba(0, 0, 0, 0.5);
          padding: 10px;
      }
      .carousel {
          max-width: 500px;
          margin: auto;
      }
      .carousel-inner img {
          max-width: 100%;
          height: auto;
      }
    </style>

    <div id="carouselExampleIndicators" class="carousel slide" data-ride="carousel">
      <ol class="carousel-indicators">
        <li data-target="#carouselExampleIndicators" data-slide-to="0" class="active"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="1"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="2"></li>
        <li data-target="#carouselExampleIndicators" data-slide-to="3"></li>
      </ol>
      <div class="carousel-inner">
        <div class="carousel-item active">
          <img class="d-block w-100" src="_static/latticefront3.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeleft3.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticeperspective3.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/latticetop3.png" alt="Top view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Top View</h5>
          </div>
        </div>
      </div>
      <a class="carousel-control-prev" href="#carouselExampleIndicators" role="button" data-slide="prev">
        <span class="carousel-control-prev-icon" aria-hidden="true"></span>
        <span class="sr-only">Previous</span>
      </a>
      <a class="carousel-control-next" href="#carouselExampleIndicators" role="button" data-slide="next">
        <span class="carousel-control-next-icon" aria-hidden="true"></span>
        <span class="sr-only">Next</span>
      </a>
    </div>