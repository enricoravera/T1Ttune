USER GUIDE
==========

The software is executed from command line, and the main entry point is the ``t1t2ne`` command. 
To get a list of available subcommands, you can run:
::

    t1t2ne --help

To get more information on a specific subcommand, you can run:
::
    
    t1t2ne <subcommand> --help

The following subcommands are available:
-   ``interactive``: an interactive command line tools that guides the user through the selection of the appropriate parameters for :math:`^15`N relaxation measurements.
-   ``ns``: a command line tool to estimate the number of scans required for a given experiment, starting from an HSQC experiment.
-   ``setuptract``: a command line tool to set up a TRACT experiment on the spectrometer, providing the appropriate delay list based on the correlation time or molecular weight.
-   ``tract``: a command line tool to fit a TRACT experiment and extract the average correlation time of the system.
-   ``makelists``: a command line tool to generate the lists for running :math:`^15`N T1 and T2 relaxation experiments given the correlation times or the molecular weight.
-   ``solventpre``: a command line tool to generate the lists for running :math:`^1`H T1 and T2 solvent PRE experiments.

The ``interactive`` subcommand
******************************

.. warning::

    The ``interactive`` subcommand requires to be run on a spectrometer workstation. The other subcommands can be run on any computer, and do not require to be run on the spectrometer.

This subcommand uses the name scheme of the two pulse sequences for T1 and T2 measurements provided with this package (hsqct1etf3gptcwg1d.rav and hsqct2etf3gptcwg1d.rav).
Assume we have a protein of 8.6 kDa (ubiquitin) and we want to measure T1 and T2 at 600 MHz.
The first experiment to be set up is the T1.
Set the sequence to hsqct1etf3gptcwg1d.rav. Calling the ``interactive`` subcommand as:
::

    t1t2ne interactive --MW 8.6 --Larmor 600

will provide the user with the maximum expected delay (``d29``), and prompt to acquire a reference experiment with a short delay (``d31 = 20 us``).
The longest delay is used for temperature compensation.
After acquiring the reference experiment, the user is prompted to provide the experiment number of the reference, phase it, and select the integration region.
The software will then provide an estimate for ``d31`` as the one that provides a reduction to about 30% of the reference intensity. 
.. note::

    This is the value that maximizes the Fisher Information criterium for a single exponential decay.

After the user has acquired the T1 experiment with the suggested delay, the software will integrate it in the same region provide a new estimate for the delay. The process is repeated until the obtained value falls within 10% of the suggested value. With the obtained T1 value, also the suggested interscan delay for the hetNOE experiment is provided.
The same process is then repeated for the T2 experiment, using the sequence hsqct1etf3gptcwg1d.rav and setting the ``l29`` and ``l31`` to the suggested initial values.
The only difference in this procedure is that the software will **NEVER** suggest to acquire a T2 experiment with a delay longer than 250 ms to avoid damage to probe, or sample, or spectrometer hardware. 

The ``ns`` subcommand
*********************

This subcommand estimates the number of scans required for a given experiment, starting from an HSQC experiment.
It is called as:
::

    t1t2ne ns --MW 8.6 --basedir <basedir> --hsqc <HSQC experiment number> --xred 10 30 70

The user needs to provide the experiment number of the HSQC, and the expected reduction in intensity for the experiment to be performed. A calculation is performed for each entry of the --xred argument.
The HSQC spectrum is phased and integrated, the number of signals under this region are counted based on the MW or provided by the user,
the SNR for each signal is then estimated and, on this basis, the number of scans required to reach a SNR of 10 at the longest delay is estimated.
The number of scans is rounded to the nearest multiple of 4 for the hetNOE experiment, and to the nearest multiple of 8 for the T1 and T2 experiments to respect phase cycling.

The ``setuptract`` subcommand
*****************************

This subcommand generates the delay list for running a TRACT experiment on the spectrometer, providing the appropriate delay list based on the correlation time or molecular weight.
It is called as:
::

    t1t2ne setuptract --MW 8.6 --Larmor 600

The user needs to provide the molecular weight or the correlation time of the system, and the Larmor frequency, unless the command is run on the spectrometer. 
The software will then generate a delay list for the TRACT experiment that covers a range of delays that are appropriate for the expected correlation time, and that maximizes the Fisher Information criterium for matching the decay rates of both :math:`\alpha` and :math:`\beta` components of the TROSY.

The ``tract`` subcommand
************************

This subcommand fits a TRACT experiment and extract the average correlation time of the system.
It is called as:
::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --integrate

.. note::

    The user needs to provide the experiment number of the TRACT experiment. If not provided, the software will look for the configuration file. If not found, it will fall back on the ``examples`` directory.
    
The software will then load the spectrum, phase it, integrate it, and fit the obtained decay curves to extract the average correlation time of the system.
Instead of integrating, the user can also select a region of the spectrum to be fitted, by providing the ``--selectregion`` argument.
The plot of the extracted correlation times as a function of the position in the spectrum is provided if the ``--plot`` argument is provided.
