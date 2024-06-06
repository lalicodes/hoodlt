.. _AnalysisExample1:

Example 1: Obtain the data from _therm and _sim
===============================================

This code reads the result of a simulation as contained in

.. code-block:: bash

    cAu1072S-Hydrocarbon-n11+cAu201S-Hydrocarbon-n11_cMgzn2_l221_a100_pDry-Ncs_uAngAmuEv_l_00p2_therm.gsd
    cAu1072S-Hydrocarbon-n11+cAu201S-Hydrocarbon-n11_cMgzn2_l221_a100_pDry-Ncs_uAngAmuEv_l_00p2_sim.gsd

The current script obtains the data in the '_therm.gsd' file and prints the average, standard deviation and
standard error. It also obtains the data in the '_sim.gsd' file


.. code-block:: python

    from hoodlt.Analysis.Collect.CollectLogData import CollectLogData
    from hoodlt.Analysis.Compute.ComputeStatistics import StatComputeTimeSeries

    file_name = 'cAu1072S-Hydrocarbon-n11+cAu201S-Hydrocarbon-n11_cMgzn2_l221_a100_ffDry-Ncs_uAngAmuEv_l_00p2'

    ff = 'dry-ncs'

    cls_log = CollectLogData([file_name], ff)
    list_quant = cls_log.quant
    print(list_quant)
    for quant in list_quant:
        mat = cls_log.get_quantity(quant)[0]
        stat = StatComputeTimeSeries(mat)
        val = stat.compute_average()
        s_err = stat.compute_std_error()
        s_dev = stat.compute_std_deviation()
        print(quant, val, s_dev, s_err)

    cls_log = CollectLogData([file_name], ff, num_frames=4, typ='sim')
    list_quant = cls_log.quant
    print(list_quant[1], cls_log.get_quantity(list_quant[1])[0,-1])
    for quant in list_quant[2:]:
        mat = cls_log.get_average(quant)[0]
        print(quant, mat[0], mat[1], mat[2])
