.. _HOODLTpredictlatticeconstant:

HOODLT prediction of lattice constants
======================================

As an example, I will illustrate the :math:`\mbox{Li}_3\mbox{Bi}` BNSL.
The experimental data is from `Talapin and Boles <https://pubs.acs.org/doi/10.1021/jacs.5b00839>`_.

This BNSL was obtained combining :math:`\mbox{Au}-\mbox{C}_{18}` with :math:`\mbox{Fe}_2\mbox{O}_3-\mbox{C}_{9}` NCs.

Here :math:`\mbox{C}_n` are thiolated (linked by a Sulphur atom) hydrocarbon chains containing :math:`n` carbons.

The diameter of :math:`\mbox{Au}` core is 4.1 nm, while the :math:`\mbox{Fe}_2\mbox{O}_3` core has 10.2 nm
The grafting densities (see TableS6 in SI) are :math:`\sigma_{Au}=5.5, \sigma_{Fe2o3}=4.5`
given in units of :math:`\mbox{nm}^{-2}`)

The **Hard Sphere** (HS) radius corresponding to both NCs is given by the OPM formula. The value of
:math:`\gamma=\frac{r_B}{r_A}` is the experimental equivalent of the same parameter considered for
packing calculations see :ref:`HOODLTpackplots`

.. math::

    r_A = 6.01 (\mbox{ nm })

    r_B = 3.41 (\mbox{ nm })

    \gamma = \frac{r_B}{r_A}=0.57 (\mbox{ dimensionless } )

where :math:`r_{A/B}` are the respectively **HS** radius,
with A= :math:`\mbox{Fe}_2\mbox{O}_3-\mbox{C}_{9}` and B= :math:`\mbox{Au}-\mbox{C}_{18}`

We now predict the packing fraction for the :math:`\mbox{Li}_3\mbox{Bi}` BNSL using the OPM and OTM
formulas it is:

.. math::

    \eta_{OPM} = \eta_{HS} = 0.55

    \eta_{OTM} = 0.84

    \eta_{experimental} = 0.88

This clearly illustrates the accuracy of the **OTM** and the breakdown of the **HS** prediction. In fact,
given the uncertainties of the experimental data, the difference of 0.04 between OTM and experiment
is within the measurement errors. With a little more work, one can predict the "fluctuation distance" between A-B NCs
within the unit cell, this gives :math:`\bar{d}_{AB}=0.98, 0.90 \mbox{ (experimental) }` in nm. The HS prediction
would be :math:`d_{AB}=2.24` and is clearly way off.

This is the script used to compute the results

.. code-block:: python

    import numpy as np
    from hoodlt.PhysChemValues.HydroCarbon import HydroCarbon
    from hoodlt.PhysChemValues.OPM import opm_val
    from hoodlt.PhysChemValues.Phys2lambxi import phys2lambxi
    import hoodlt.OTM.Li3Bi_OTM as OTMLat

    # diameters of the cores
    diam_fe2o3 = 10.2
    sigma_fe2o3 = 4.5
    diam_au = 4.1
    sigma_au = 5.5
    # experimental pf
    pf_experimental = 0.88

    # relevant hydrocarbon objects
    hc18 = HydroCarbon(18)
    hc9 = HydroCarbon(9)

    # compute radius
    rad_au = diam_au*0.5
    rad_fe2o3 = diam_fe2o3*0.5

    # convert the parameters into the dimensionless quantities lambda, xi
    auc18 = phys2lambxi(hc18.max_length(), hc18.area_dens_max, rad_au, sigma_au)
    fe2o3c9 = phys2lambxi(hc9.max_length(), hc9.area_dens_max, rad_fe2o3, sigma_fe2o3)
    opm_fe2o3c9
    # get the opm = hard sphere radius
    opm_auc18 = rad_au*opm_val(auc18[0], auc18[1])
    opm_fe2o3c9 = rad_fe2o3*opm_val(fe2o3c9[0], fe2o3c9[1])

    # the value gamma is just the ratio of the two hard sphere radius
    gam = opm_auc18/opm_fe2o3c9

    print('particle A radius', opm_fe2o3c9)
    print('particle B radius', opm_auc18)
    print('gamma ', gam)
    print('\n')

    # we choose a [2,2,2] lattice with lattice constant = 1, the lattice constant here is in arbitrary units
    # the value of gamma needs to be the one predicted above
    lat = OTMLat.OTMLatLi3BiBase16([2,2,2], 1, gam)

    # those are the packing fraction predictions
    print('opm pf', lat.pf())
    print('otm pf', lat.otm_observables().pf)
    print('exp   ', pf_experimental)

