.. _NcsExample5:

Example 5: build a simulation of water with the spce model
========================================================

Here is the script to model water within the spce model, which is defined in the opls-aa force field.

.. code-block:: python

    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelsolvent.SpceWater import SpceWater

    # this is the name as the forcefield you have chosen to use to build your configuration
    ff  = 'opls-aa'

    # number of water molecules
    num_w = 1000
    # grid is in units of nanometers
    h_grid = 0.3
    L_max = 3
    box = [L_max, L_max, L_max, 0.0, 0.0, 0.0]

    # this is the configuration builder object
    builder = ConfigurationBuilder()

    # add solvent
    builder.add_solvent([SpceWater(ff)]*num_w, [h_grid, h_grid, h_grid] , box)

    # set the alias for the configuration
    builder.set_alias('spce'+'_'+str(num_w))

    # get the configuration we just built
    conf = builder.get_configuration()

    # save the configuration, for simulating later
    save_config(conf)
**Expected Simulation Views**


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
          <img class="d-block w-100" src="_static/buildfront2.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/buildleft2.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/buildperspective2.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/buildtop2.png" alt="Top view">
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