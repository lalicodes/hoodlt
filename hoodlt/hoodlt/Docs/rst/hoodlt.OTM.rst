The Orbifold Topological Model (OTM)
====================================

The OTM model computes the lattice structure by describing a nanocrystal as a soft Skyrmion that supports texture defects. The key elements are the values of :math:`\bar{\gamma}`, the nanocrystal diameters :math:`{\bar r}_A, {\bar r}_B` and the OTM packing fraction :math:`\bar{\eta}`. The model is fully described in two references:

`Topological structure prediction in binary nanoparticle superlattices, Travesset, A.;Soft Matter, 2017,13, pp 147-157 <http://pubs.rsc.org/en/content/articlelanding/2017/sm/c6sm00713a#!divAbstract>`_

`Soft Skyrmions, Spontaneous Valence and Selection Rules in Nanoparticle Superlattices, Travesset, A.; ACS Nano, 2017, 11 (6), pp 5375–5382 <http://pubs.acs.org/doi/abs/10.1021/acsnano.7b02219>`_

Please cite one of the two papers(or both) if you use OTM.

Examples
--------

Here it is how it is used. In order to create
a :math:`\mbox{MgZn}_{2}` OTM lattice object, we use the following
snippet.

.. code-block:: python

   import hoodlt.OTM.MgZn2_OTM as Om

        l_box = 3
        # define a box with 3 x 3 x 3 unit cells

        a_nn = 1
        # minimum separation between A-particles 

        gamma = 0.8
        # ratio of the diameter of B and A particle

        mgz = Om.OTMLatMgZn2Base12(l_box, a_nn, gamma)
        # create lattice object
        
        print(mgz.pf()) 
        # hard sphere packing fraction 

        obs = mgz.observables()
        # obtain the OTM object

        print(obs.ind)
        # print the OTM index

        print(obs.pf)
        # OTM packing fraction

        print(obs.ratio_a)
        # ratio of a
        
        print(obs.gamma_bar)
        # value of gamma_bar

OTM Classes
-----------
        
.. toctree::

    hoodlt.OTM_list.rst
