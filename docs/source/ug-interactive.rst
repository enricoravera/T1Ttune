.. _sec;interactive:

The ``interactive`` subcommand
******************************

.. warning::

    The ``interactive`` subcommand requires to be run on a spectrometer workstation.
    The other subcommands can be run on any computer, included personal and processing stations.

This subcommand uses the name scheme of the two pulse sequences for T\ :sub:`1`\  and T\ :sub:`2`\  measurements provided with this package (`hsqct1etf3gptcwg1d.rav` and `hsqct2etf3gptcwg1d.rav`).

Let us assume we have a protein of 8.6 kDa (ubiquitin) and we want to measure T\ :sub:`1`\  and T\ :sub:`2`\  at 600 MHz.


-   The first experiment to be set up is the T\ :sub:`1`\ .  
    Set the sequence to `hsqct1etf3gptcwg1d.rav`. 
-   Calling the ``interactive`` subcommand as:
    ::

        t1t2ne interactive --MW 8.6 --Larmor 600

    will provide the user with the maximum expected delay (``d29``), and prompt to acquire a reference experiment with a short delay (``d31 = 20 us``).
    The longest delay is used for temperature compensation.

-   After acquiring the reference experiment, the user is prompted to provide the experiment number of the reference, phase it, and select the integration region.
    The instructions on how to use the GUIs can be found in :func:`klassez.processing.interactive_phase_1D` and :func:`klassez.anal.integrate` respectively.
   
    The software will then provide an estimate for ``d31`` as the one that provides a reduction to about 30% of the reference intensity. 

    .. note::

        This is the value that maximizes the `Fisher Information criterium`_ for a single exponential decay.

    .. _Fisher Information criterium: https://en.wikipedia.org/wiki/Fisher_information

-   After the user has acquired the T\ :sub:`1`\  experiment with the suggested delay, the software will integrate it in the same region provide a new estimate for the delay.
    The process is repeated until the obtained value falls within 10% of the suggested value. 
  
With the obtained T\ :sub:`1`\  value, also the suggested interscan delay for the hetNOE experiment is provided.

The same process is then repeated for the T\ :sub:`2`\  experiment, using the sequence `hsqct2etf3gptcwg1d.rav` and setting the ``l29`` and ``l31`` to the suggested initial values.
The only difference in this procedure is that the software will **NEVER** suggest to acquire a T\ :sub:`2`\  experiment with a delay longer than 250 ms to avoid damage to probe, or sample, or spectrometer hardware. 

The software accepts several command line options to customize the procedure:

-   Parameters related to the acquisition of the experiments:
    -   ``--basedir``: the base directory of the experiments. This overrides the default current directory. If not set, the user will be prompted to enter it at the beginning of the procedure. The software will check that the specified directory exists.
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

.. seealso:: 

   :class:`t1t2ne.scripts.t1t2ne_interactive`
