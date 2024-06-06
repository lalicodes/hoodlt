.. _AnalysisExample5:

Example 5: Cluster Potential of Mean Force
==========================================

The following script can be used with the output of a simulation similar to :ref:`SimulationExample4` done on a
configuration of nanoparticles similar to :ref:`NcsExample3` to compute the potential of mean force for the entire
system, and compare that with the potential associated with the sum of pair potentials in the system.


.. code-block:: python

    from hoodlt.Analysis.Plot.PlotClusterPMF import PlotClusterPMF
    from hoodlt.Analysis.Compute.ComputeClusterPMF import ComputeClusterPMF
    from hoodlt.Data.Modelconfigurations.CommonConfigurations.PlanarConfigData import PlanarConfigData

    ff = "dry-ncs"

    rvals = [59, 57, 55, 51, 47, 43, 39, 35, 31, 30, 29, 28, 27, 26, 25]

    file_beg = "Au201-Hydrocarbon-n11_cP2_r"
    file_end = "pNcs-In-Solvent_uAngAmuEv"
    files = []
    for i in range(len(rvals)):
	    files.append(file_beg + str(rvals[i]) + file_end)

    # this is the configuration info object used to make the configuration
    config = PlanarConfigData(2, 55)

    # create compute and plot objects
    comp = ComputeClusterPMF(files, rvals, config, ff, 'CTR-CTR1', 'CTR-CTR2', pair_pmf_file='201-12_potential.txt')
    plotter = PlotClusterPMF(comp)

    # these functions plot the results of the calculation
    # even if the calculation is not done first explicitly by the compute class, the plot class will do the calculation internally
    plotter.plot_pmf()
    plotter.write_pmf_plot_data()

    plotter.plot_many_body_effects(0, 55, 201, 20.3, 80, 12)
    plotter.write_many_body_plot_data()

    plotter.plot_pair_pmf_comparison()
    plotter.write_pair_pmf_comparison_plot_data()
