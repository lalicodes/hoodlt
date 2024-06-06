.. _NcsExample8:

Example 8: Perovskite Nanocube Pair, Planar Square Lattice
==========================================================

Here we construct Perovskite nanocube pair systems and planar simple square lattice systems.

1. Building a nanocube pair is similar to building an NC pair as in :ref:`NcsExample2`.
2. Building a planar square lattice is similar to building a lattice as in :ref:`NCsExample4`.

.. code-block:: python

    import numpy as np
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelnanoparticles.Perovskite import Perovskite
    from hoodlt.Data.Modelligands.OleateLigand import OleateLigand
    from hoodlt.Data.Modelligands.OleylammoniumLigand import OleylammoniumLigand
    from hoodlt.Lattices.Lat_dim2_1.ss_2d_lat_Mixt import LatSquare as square
    from hoodlt.Data.ProcessConfigurations.Squeeze import Squeeze
    from hoodlt.Data.Modelconfigurations.Saver import save_config

    # Initialize configuration and forcefield
    forcefield = 'opls_dry-ncs_mix'
    l_size = 2
    builder = ConfigurationBuilder()

    # Create nanocube pair
     core = Perovskite(forcefield=forcefield, l_size=l_size)
    lig1 = OleateLigand(forcefield, straight=True)
    lig2 = OleylammoniumLigand(forcefield, straight=True)
    ligands = [lig1]*core.graft_num[0] + [lig2]*core.graft_num[1]
    builder.add_nc(core, ligands, [-35, 0, 0])
    builder.add_nc(core, ligands, [35, 0, 0])
    builder.add_bond(0, 1, 'CTR-CTR1')

    # Configure and build planar square lattice
    lat = square(l_value=3, a_nn_e=50)
    lat.resize_box(3*50)
    builder.add_nc(core, ligands, [0, 0, 0])
    nc = builder.conf.particles[0]
    list_nc = Squeeze(forcefield, list_nc=[nc], list_radius=[23]).squeeze()
    builder.build_lattice_from_reinit(list_nc, lat)

    # Set configuration parameters and save
    builder.set_box(300)
    builder.set_alias('Pair & Lattice Simulation')
    save_config(builder.conf, filename='simulation.gsd')






**Expected Simulation Views**

**First One:**

.. raw:: html

    <link rel="stylesheet" href=""https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
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
          <img class="d-block w-100" src="_static/peronanofront.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanoleft.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanoperspective.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanotop.png" alt="Top view">
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


**Second One:**

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
          <img class="d-block w-100" src="_static/peronanofront2.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanoleft2.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanoperspective2.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/peronanotop2.png" alt="Top view">
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
