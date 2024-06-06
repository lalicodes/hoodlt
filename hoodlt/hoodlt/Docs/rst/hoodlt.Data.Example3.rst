.. _NcsExample3:

Example 3: NC Cluster
=====================

Examples 1 and 2 were simple systems that only had 1 or 2 nanoparticles; it was easy to just quickly add them to the
system without much thought as to what their positions would be and what indexes would be bonded together. For larger
systems of nanoparticles, this is not the case. That is why HOODLT uses a set of objects located in
:mod:`hoodlt.Data.ModelConfigurations.CommonConfigurations` to store lists of positions and bond indexes to quickly
and easily build more complex configurations of nanoparticles that we would often like to build. An example using one
of is shown below, where we build an icosahedron consisting of a central nanoparticle and twelve
nearest neighbors.

.. code-block:: python

    from hoodlt.Data.Modelnanoparticles.TO201 import TO201
    from hoodlt.Data.Modelnanoparticles.TO140 import TO140
    from hoodlt.Data.Modelconfigurations.Saver import save_config
    from hoodlt.Data.Modelconfigurations.ConfigurationBuilder import ConfigurationBuilder
    from hoodlt.Data.Modelligands.HydrocarbonLigand import HydrocarbonLigand as Hydrocarbon
    from hoodlt.Data.Modelconfigurations.CommonConfigurations.IcosConfigData import IcosConfigData

    forcefield = "ncs-in-solvent"

    core = TO201(forcefield)
    core2 = TO140(forcefield)
    lig = Hydrocarbon(repeats=11, ff=forcefield)
    lig2 = Hydrocarbon(repeats=8, ff=forcefield)

    # this is a configuration data object, which makes building larger configurations much simpler
    conf_data = IcosConfigData(dist_from_origin=70, nc_at_origin=True)

    builder = ConfigurationBuilder()

    # here we loop overif True: the list of positions held by the configuration data object, to quickly add the ncs
    # here we are adding a single Au201 11 repeat Hydrocarbon NC at the origin and many Au140 8 repeat Hydrocarbon NCs
    # at the surrounding positions
    for i, pos in enumerate(conf_data.nc_positions):
            if i == 0:
                    builder.add_nc(core, [lig]*core.graft_num, pos)
            else:
                    builder.add_nc(core2, [lig2]*core2.graft_num, pos)

    # here we loop over the list of bond indexes held by the configuration data object, to quickly add the bonds
    # here we are adding all bonds connected to the Au201S NC as CTR-CTR1 bonds, and all the other bonds are CTR-CTR2
    ctr_bond_list = ['CTR-CTR1', 'CTR-CTR2']
    for i, j, k in conf_data.bond_list:
            builder.add_bond(i, j, ctr_bond_list[k])


    # you need to define a box containing the nc
    builder.set_box(300)

    # set the name
    builder.set_alias("Icosahedron")

    conf = builder.get_configuration()

    save_config(conf)


The :mod:`hoodlt.Data.Modelconfigurations.CommonConfigurations` directory includes classes for building many types of
configurations, listed here: :ref:`NanocrystalConfigs` The meaning of the output files from this script are similar to
that of :ref:`NcsExample1`.

Note that we had to define the box before saving the configuration.

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
          <img class="d-block w-100" src="_static/ncclusterfront.png" alt="Front view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Front View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/ncclusterleft.png" alt="Left view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Left View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/ncclusterperspective.png" alt="Perspective view">
          <div class="carousel-caption d-none d-md-block">
            <h5>Perspective View</h5>
          </div>
        </div>
        <div class="carousel-item">
          <img class="d-block w-100" src="_static/ncclustertop.png" alt="Top view">
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