.. _NcsExample2:

Example 2: NC pair
==================

Example 1 was a simple system that only had a single nanoparticle; but we can do a lot more with the
:mod:`ConfigurationBuilder` object

.. code-block:: python

    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand as Hydrocarbon
    from hoodlt.Data.Modelsolvent.Toluene import Toluene

    # we will add num_solvent number of solvent molecules (toluene in this example)
    num_solvent = 1000
    # each solvent will be placed in a cubic box of h_grid size
    h_grid=10
    # separation between the two nanoparticles
    l_dist = 55.0
    # simulation box
    box = [3*l_dist, 3*l_dist, 3*l_dist, 0.0, 0.0, 0.0]

    # forcefield
    forcefield = "ncs-in-solvent"

    # core
    core = TO201(forcefield)
    # ligand
    lig = Hydrocarbon(repeats=11, ff=forcefield)

    # this is the solvent, which will be added after we add all the nanoparticles
    solv = Toluene(forcefield)

    builder = ConfigurationBuilder()

    # here we add multiple nanoparticles to our configuration
    ind1 = builder.add_nc(core, [lig]*core.graft_num, [-0.5*l_dist, 0, 0])
    ind2 = builder.add_nc(core, [lig]*core.graft_num, [0.5*l_dist, 0, 0])

    # we can also add bonds between bondable entities in our configuration
    builder.add_bond(ind1, ind2, 'CTR-CTR1')

    # Now that all the nanoparticles have been added to the system, we add the solvent
    builder.add_solvent(num_solvent*[solv], [h_grid, h_grid, h_grid], box)

    builder.set_alias("Pair")

    conf = builder.get_configuration()

    save_config(conf)

The meaning of the output files from this script are similar to that of :ref:`NcsExample1`.


A Note on Solvent Addition
--------------------------

The h_grid variable define the dimensions of the cell that each solvent will be contained
in. When you choose to add the solvent, you will probably already know how many solvents you want to add and at what
density you want to add them. Getting the tolerance right is
a matter of trial and error, and you may need to try a few values before you get it to fill the box in the way you want.

If upon trying to save the configuration after you added the solvent, and HOOMD complains about particles which are
outside the box, you will need to increase the grid size.

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
          <img class="d-block w-100" src="_static/ncfront.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/ncleft.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/ncperspective.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/nctop.png" alt="Top view">
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

    <script>
      // Ensure the carousel does not auto-cycle
      document.addEventListener('DOMContentLoaded', function() {
        $('.carousel').carousel({
          interval: false
        });
      });
    </script>
