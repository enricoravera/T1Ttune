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

.. seealso:: 

   :class:`t1t2ne.scripts.t1t2ne_interactive`
