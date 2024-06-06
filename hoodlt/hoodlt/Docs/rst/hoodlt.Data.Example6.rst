.. _NcsExample6:

Example 6: Substrates
=====================


HOODLT also has support for adding substrates into the simulation box. Be careful and make sure that the length of your
substrate is not too much bigger than 1/3 the length of the simulation box, or else HOOMD will fail to start the
simulation. To avoid this error, we set the box explicitly in the script.

.. code-block:: python

    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelsubstrates.SquareSubstrate import SquareSubstrate
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelsolvent.Toluene import Toluene
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand as Hydrocarbon
    from hoodlt.Data.Modelconfigurations.Saver import save_config

    forcefield = "ncs-in-solvent"

    core = TO201(forcefield)
    lig = Hydrocarbon(11, forcefield)

    # here we build the substrate object
    subs = SquareSubstrate(forcefield, grid_dim=[20, 20], spacing=1.2, particle_type='AuSub')

    builder = ConfigurationBuilder()

    # this is how we add substrates to a configuration
    builder.add_substrate(subs, [0, 0, 0])

    # you should only set the box explicitly if you are working with substrates. Otherwise, let the ConfigurationBuilder
    # do it for you
    builder.set_box([111*1.2*3]*3)

    # add the ncs
    builder.add_nc(core, [lig]*core.graft_num, [-27.5, 0, 27.5])
    builder.add_nc(core, [lig]*core.graft_num, [27.5, 0, 27.5])

    conf = builder.get_configuration()
    save_config(conf)

Here you should see the files:

.. code-block:: bash

    Au201-Hydrocarbon-n11_bAuSub20-20_ffNcs-In-Solvent_uAngAmuEv_bonds.json
    Au201-Hydrocarbon-n11_bAuSub20-20_ffNcs-In-Solvent_uAngAmuEv.pickle
    Au201-Hydrocarbon-n11_bAuSub20-20_ffNcs-In-Solvent_uAngAmuEv_restart.gsd

**Expected simulation views:**

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
          <img class="d-block w-100" src="_static/substratefront.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/substrateleft.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/substrateperspective.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/substratetop.png" alt="Top view">
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
