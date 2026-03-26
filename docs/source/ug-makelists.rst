.. _sec;makelists:

The ``makelists`` subcommand
****************************

This subcommand generates the lists for running :sup:`15`\ N T\ :sub:`1`\  and T\ :sub:`2`\  relaxation experiments given the correlation times or the molecular weight.
It is called as:

::
    
    t1t2ne makelists --MW 8.6 --Larmor 600


or, following the ubiquitin TRACT experiment described below in :ref:`sec;tract`:

::

    t1t2ne makelists --tau 4.57e+00 --S2 0.90 --Larmor 600


Likewise, after the synuclein TRACT experiment described in :ref:`sec;tract`:

::

    t1t2ne makelists --tau 8.80e+00 1.49e+00 --S2 0.15 0.38 --Larmor 600


The software accepts several command line options to customize the procedure:

-   Parameters related to the acquisition of the experiments:
    -   ``--xred``: the percent residual signal for the longest delay of the experiments. If not set, it will defaulto to 5% and 10% for T\ :sub:`1`\ and T\ :sub:`2`\  respectively.
    -   ``--nT``: the number of increments in the suggested vdlist for the T\ :sub:`1`\  experiment and T\ :sub:`1`\  T\ :sub:`2`\  experiments. Default is 8. If more than one value is provided, the first one will be used for the T\ :sub:`1`\  experiment, and the second one for the T\ :sub:`2`\  experiment.
    -   ``--logscale``: whether to use a logarithmic scale for the vdlist and vclist. Default is False, which means a linear scale will be used.
    -   ``--randomize``: whether to randomize the order of the values in the lists. Default is False.
    -   ``--Larmor``: the Larmor frequency of the spectrometer, in MHz. It will anyway be read from the configuration file of the spectrometer.
    -   ``--B0``: the magnetic field strength of the spectrometer, in Tesla. It will anyway be read from the configuration file of the spectrometer.
    -   ``--large``: whether to create the lists for the "large" sequence, which is optimized for short T2 times. Default is False.
    -   ``--small``: whether to create the lists for the "small" sequence, which is optimized for long T2 times. Default is False.
-   Parameters related to the estimation of the correlation times and order parameters:
    -   ``--S2``: the Lipari-Szabo order parameter S2 to use for the calculation of tau_c. In IDP mode, two values should be provide, else only one
    -   ``--MW``: the molecular weight of the protein in kDa, to be used for estimating the longest correlation time (tau_slow in the IDP model).
    -   ``--tau``: the correlation times to use for the estimation of the initial delays. In IDP mode, two values should be provided, else only one. If not set, they will be estimated from the MW. Otherwise, default T\ :sub:`1`\ and T\ :sub:`2`\  values will be used for the non-IDP model, and the three-tau model described in `Bolognesi et al.`_ for the IDP model.
    -   ``--idp``: whether the system under study is an intrinsically disordered protein (IDP). If set, the software will use the three-tau model described above. Default is False.
        -   If set, it will also select a longer delay in the cpmg for the T\ :sub:`2`\  experiment, optimized for the expected relaxation behavior of IDPs as described by `Bolognesi et al.`_.
        -   If set, it will also select a longer delay for the T\ :sub:`1`\  experiment, to avoid excessive heating, because cross-correlations are less important in IDPs.
        -   the ``--MW`` argument will be used to compute the slowest correlation time.
        -   ``--corr_window_idp``: the length of the correlation window for the IDP model, in residues. Default is 20.
    -   ``--T``: the temperature in Kelvin, to be used for estimating the taus from MW. Default is 298 K.
-   Parameters related to relaxation:
    -   ``--r``: the length of the 1H-15N bond in Angstroms. Default is 1.02 A.
    -   ``--Deltasigma``: the chemical shift anisotropy of the 15N nucleus in ppm. Default is -160 ppm.
    -   ``--nucs``: the nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.

.. Bolognesi et al: https://doi.org/10.1016/j.pnmrs.2025.101577

.. rubric:: TODO

Chiamata generale con le opzioni 

Sezione esempi

Chiamate varie in SEEALSO
