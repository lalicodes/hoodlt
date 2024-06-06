.. _Analysis:

Analysis of HOODLT Output Files
===============================

Introduction
------------

The purpose of these modules is to analyze the output data from any HOODLT :ref:`HOODLTSimulations`. This process
can always be broken down into 3 steps:

1. **Data Collection (Collect):** Gather raw output data and transform it into actionable insights.
2. **Data Processing (Compute):** Unleash the computational power to crunch numbers and extract valuable information.
3. **Object Analysis (Analyze):** Dive deep into individual simulation objects to uncover intricate details.
4. **Visualization (Plot):** Craft visually stunning representations of your analysis results.
There is a submodule for each of these categories, described below

.. toctree::
   :maxdepth: 2
   :caption: Explore Submodules

   hoodlt.Analysis.subpkgs

Usage
-----

The analysis framework provides various levels of support. While predefined plotting and file writing formats are available in :mod:`hoodlt.Analysis.Plot`, plotting data directly from objects in :mod:`hoodlt.Analysis.Compute` may sometimes be more straightforward. Additionally, you can inspect simulation output data before performing rigorous calculations by utilizing objects in :mod:`hoodlt.Analysis.Collect`.

.. toctree::
   :maxdepth: 1
   :caption: Common Usage Examples

   hoodlt.Analysis.Example1
   hoodlt.Analysis.Example2
   hoodlt.Analysis.Example3
   hoodlt.Analysis.Example4
   hoodlt.Analysis.Example5
