.. _AnalysisExample2:

Example 2: Average Nanoparticles
================================

The following script can be used to compute average positions of ligands within nanoparticles over a
range of frames from a gsd. Averaging nanoparticles allows to analyze vortex textures.
Additional analysis is needed to further determine which ligands belong to the vortices.

.. code-block:: python

    from hoodlt.Analysis.Collect.Trajectory import Trajectory
    from hoodlt.Analysis.Compute.ComputeAverageNCs import ComputeAverageNCs

    f_name = "Au201-Hydrocarbon-n11_cSingle_pNcs-In-Solvent_uAngAmuEv_restart"

    trj = Trajectory(f_name, f_name)

    comp = ComputeAverageNCs(trj)

    # computes the average of the first nanoparticle over the given number of frames
    nc = comp.compute_average_particle(0)
    conf = nc.configuration()
    conf.write_gsd("nc0_average_220-240")

After running this script, you will find a file

.. code-block:: bash

    nc0_average_220-240.gsd

that contains the nanoparticle whose ligand positions are averaged over the entire run.