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

-   ``interactive``: an interactive command line tools that guides the user through the selection of the appropriate parameters for :math:`^{15}` N relaxation measurements.
-   ``ns``: a command line tool to estimate the number of scans required for a given experiment, starting from an HSQC experiment.
-   ``setuptract``: a command line tool to set up a TRACT experiment on the spectrometer, providing the appropriate delay list based on the correlation time or molecular weight.
-   ``tract``: a command line tool to fit a TRACT experiment and extract the average correlation time of the system.
-   ``makelists``: a command line tool to generate the lists for running :math:`^{15}` N T1 and T2 relaxation experiments given the correlation times or the molecular weight.
-   ``solventpre``: a command line tool to generate the lists for running :math:`^{1}` H T1 and T2 solvent PRE experiments.

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

This subcommand fits a TRACT experiment and extracts the average correlation time of the system. It uses the algebraic analysis described in Robson et al. (2021), `doi:10.1007/s10858-021-00379-5`_.
    
    .. _doi:10.1007/s10858-021-00379-5: https://doi.org/10.1007/s10858-021-00379-5

It is called as:
::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --integrate

.. note::

    The user needs to provide the experiment number of the TRACT experiment. If not provided, the software will look for the configuration file. If not found, it will fall back on the ``examples`` directory.
    
The software will then load the spectrum, phase it, integrate it, and fit the obtained decay curves to extract the average correlation time of the system.
Instead of integrating, the user can also select a region of the spectrum to be fitted, by providing the ``--selectregion`` argument.
The plot of the extracted correlation times as a function of the position in the spectrum is provided if the ``--plot`` argument is provided.

.. figure:: _static/tract_ubiquitin.png
    :name: tract_ubiquitin
    :width: 80.0%

    The result of the TRACT analysis for ubiquitin at 600 MHz in the ``--selectregion`` mode. The raw data are provided in the ``examples`` directory.

If the ``--idp`` option is used, the software will **NOT** compute the correlation time but the order parameter :math:`S^2_{int}` of the intermediate motion. 
The underlying assumption is that the global reorientation of an IDP corresponds to the global reorientation of a globular protein of the same weight (hence the need for the ``--MW`` argument), but has a low order parameter, and that a second motion, corresponding to 20-30 residues (intermediate motion) is present. This model is described in the work of Razaei-Ghaleh et al. (2018), `doi:10.1002/anie.201808172`_.
    
    .. _doi:10.1002/anie.201808172: https://doi.org/10.1002/anie.201808172 

In this case the call is, for instance for synuclein at 600 MHz:
::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --idp --MW 14.4 --selectregion --plot

The results are shown in the following figure.

.. figure:: _static/tract_synuclein.png
    :name: tract_synuclein
    :width: 80.0%

    The result of the TRACT analysis for the IDP synuclein at 600 MHz in the ``--selectregion`` mode. The raw data are provided in the ``examples`` directory.

At the end of the analysis, the software provides a command to generate the lists for running T1 and T2 experiments based on the obtained correlation time.

The ``makelists`` subcommand
****************************

This subcommand generates the lists for running :math:`^15`N T1 and T2 relaxation experiments given the correlation times or the molecular weight.
It is called as:

::
    
    t1t2ne makelists --MW 8.6 --Larmor 600

or, following the ubiquitin TRACT experiment:

::
    t1t2ne makelists --tau 4.57e+00 --S2 0.90 --Larmor 600


Likewise, for the synuclein TRACT experiment:

::

    t1t2ne makelists --tau 8.80e+00 1.49e+00 --S2 0.15 0.38 --Larmor 600

The user needs to provide the molecular weight or the correlation times of the system, and the Larmor frequency, unless the command is run on the spectrometer.
By default, lists are generated with linear spacing, however, logarithmic spacing is more appropriate and can be activated with the ``--logspace`` argument.

The ``solventpre`` subcommand
*****************************

This subcommand generates the lists for running :math:`^1` H T1 and T2 solvent PRE experiments. The expected relaxation rates are computed with the Freed Outer Sphere model. 
In this model, the spectral density function is computed as:

.. math::

    \tau_D = \frac{d^2}{D_{target} + D_{cosolute}} \breqn
    z = \sqrt{2 |\omega| \tau_D + \frac{\tau_D}{\tau_1}} \breqn
    J_{outer}(\omega) = \frac{2}{5} \frac{1 + \frac{5 z}{8} + \frac{z^2}{8}}{1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{4 z^4}{81} + \frac{z^5}{81} + \frac{z^6}{648}}
        
where: 
    -   :math:`d` is the distance of closest approach of the paramagnetic center and the nucleus in meters
    -   :math:`D_{target}` is the diffusion coefficient of the target molecule in m^2/s
    -   :math:`D_{cosolute}` is the diffusion coefficient of the cosolute molecule in m^2/s
    -   :math:`\tau_1` is the correlation time for electron relaxation in seconds
    -   :math:`\omega` is the angular frequency where the spectral density is evaluated, in rad/s

Equations 6.42, 6.48 and 6.50 in `Bertini et al. 2016`_.
 
    .. _Bertini et al. 2016: https://www.sciencedirect.com/science/chapter/monograph/pii/B9780444634368000065

This spectral density function is then used to compute the expected PRE rates, with optional transient zero-field splitting (ZFS) contribution to the electron relaxation.
    Default values are for 1 mM Gd-DOTA and a protein of 10 kDa at room temperature. The default values correspond to Gd-DOTA, and are taken from `Li et al. 2002`_.
    
    .. _Li et al. 2002: https://pubs.acs.org/doi/full/10.1021/ic0200390

.. math:: 

    k_{outer} = \frac{16\pi}{81} \left( \frac{\mu_0}{4\pi} \right)^2 h^2 \gamma_n^2 \gamma_e^2 S(S+1) f N_A \frac{c}{d(D_{target}+D_{cosolute})}\breqn
    R_1 = k_{outer} [7 J_{outer} (\omega_S) + 3 J_{outer} (\omega_I)]\breqn
    R_2 = \frac{1}{2} k_{outer} [13 J_{outer} (\omega_S) + 3 J_{outer} (\omega_I) + 4 J_{outer} (0)]

Where:
    -   :math:`\mu_0` is the vacuum permeability in T m/A
    -   :math:`h` is the Planck constant in J s
    -   :math:`\gamma_n` is the gyromagnetic ratio of the nucleus in rad/s/T
    -   :math:`\gamma_e` is the gyromagnetic ratio of the electron in rad/s/T
    -   :math:`S` is the electron spin quantum number
    -   :math:`N_A` is the Avogadro's number in mol^-1
    -   :math:`c` is the concentration of the paramagnetic cosolute in mol/m^3
    -   :math:`f` is the fraction of the sphere around the nucleus that is accessible to the paramagnetic cosolute.
    
Optionally, the transient zero-field splitting (ZFS) contribution to electron relaxation is included with the Bloembergen-Morgan model:

.. math::

    R_{1e} = \frac{2}{15} \Delta_t^2 \left(4S(S+1) - 3\right) \left(J(\omega_e) + 2J(2\omega_e)\right)\breqn
    R_{2e} = \frac{1}{15} \Delta_t^2 \left(4S(S+1) - 3\right) \left(3J(0) + 5J(\omega_e) + J(2\omega_e)\right)

Where:
    -   :math:`\Delta_t` is the transient zero-field splitting, provided in inverse centimeters and translated automatically to Hz
    -   :math:`\omega_e` is the electron Larmor frequency in rad/s

It can be called providing only the field as:
::

    t1t2ne solventpre --Larmor 600

Alternatively, the user can provide the expected intrinsic T1 of the system and/or the expected linewidth:
::

    t1t2ne solventpre --Larmor 600 --T1 1 --lw 10

The expected linewidth can also be estimated from the molecular weight of the system:
::

    t1t2ne solventpre --Larmor 600 --MW 8.6

By default, the software will suggest a list of 8 delays for T1 and 2 points for T2, as described by Iwahara, Tang, and Clore in `10.1016/j.jmr.2006.10.003`_. The number of points can be changed with the ``--nT`` argument.
    
    .. _10.1016/j.jmr.2006.10.003: https://doi.org/10.1016/j.jmr.2006.10.003


**Important**: The delays are calculated under the assumption that the pulse sequence used for T2 is the one found in figure 1 of `10.1016/j.jmr.2006.10.003`_.:




