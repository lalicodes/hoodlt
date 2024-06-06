Explicit Potentials 
===================

Examples
--------

Here it is how it is used. In order to create
a potential object, we would use the following
snippet.

.. code-block:: python

        import hoodlt.Potentials.general_LJ as Lj

        alat = 1
        # length units
        cut_off = 4.1
        # cut-off values
        p_exp = 12
        # p- exponent
        q_exp = 6
        # q- exponent

        vlj = Lj.LJGen(p_exp, q_exp, 1.0, 1.0, cut_off)
        # create the potential object, normaly used in conjunction with other hoodlt functions

Generalized Lennard Jones 
-------------------------

.. automodule:: hoodlt.Potentials.general_LJ
    :members:
    :undoc-members:
    :show-inheritance:

Inverse power law with cutoff
-----------------------------

.. automodule:: hoodlt.Potentials.overrp_pot_cutoff
    :members:
    :undoc-members:
    :show-inheritance:

Mixture of inverse power laws
-----------------------------

.. automodule:: hoodlt.Potentials.overrp_pot_cutoff_Mixt
    :members:
    :undoc-members:
    :show-inheritance:

Mixture of power laws and nulls
-------------------------------

.. automodule:: hoodlt.Potentials.overrp_pot_cutoff_vanish_Mixt
    :members:
    :undoc-members:
    :show-inheritance:


