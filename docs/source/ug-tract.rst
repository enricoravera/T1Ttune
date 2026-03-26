.. _sec;tract:

The ``tract`` subcommand
************************

This subcommand fits a TRACT experiment and extracts the average correlation time of the system. 
It uses the algebraic analysis described in `Robson et al. (2021)`_, and that is available on GitHub as `nomadiq/TRACT_analysis`_.
    
.. _Robson et al. (2021): https://doi.org/10.1007/s10858-021-00379-5
.. _nomadiq/TRACT_analysis: https://github.com/nomadiq/TRACT_analysis/tree/master

The script is invoked as:
::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --integrate

The TRACT experiment to analyzed is interpreted to be in the `<basedir>/<TRACT experiment number>` directory.

.. note::

    The user `must` provide the experiment number of the TRACT experiment. If not provided, the software will look for the configuration file. 
    If it is not found, it will fall back on the `examples` directory.
    
The algorithm works as follows.

-    The software will load the spectrum and process it.
-    phase it, integrate it, **VEDI TODO**
-    The decay curves obtained from the integrated regions are then fitted to extract the average correlation time of the system.
-    The values can be visualized as function of the position in the spectrum by providing the ``--plot`` argument.

.. rubric:: TODO

Instead of integrating, the user can also select a region of the spectrum to be fitted, by providing the ``--selectregion`` argument, and the plot of the extracted correlation times as a function of the position in the spectrum is provided if the ``--plot`` argument is provided. **QUESTA VA SPIEGATA MEGLIO**


::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --selectregion --plot

---

The results are shown in the following figure.

.. figure:: _static/tract_ubiquitin.png
    :name: tract_ubiquitin
    :width: 80.0%

    The result of the TRACT analysis for ubiquitin at 600 MHz in the ``--selectregion`` mode. The raw data are provided in the ``examples`` directory.

**TODO** `questo è importante e va detto prima`
At variance with the original paper (`Robson et al. (2021)`_), we foresee the use for Intrinsically Disordered Proteins (IDPs). 
The underlying assumption is that the global reorientation of an IDP corresponds to the global reorientation of a globular protein of the same weight (hence the need for the ``--MW`` argument), but has a low order parameter, and that a second motion, corresponding to 10-20 residues (intermediate motion) is present. This model is described in the work of `Razaei-Ghaleh et al. (2018)`_.
Therefore, when the ``--idp`` option is used, the software will **NOT** compute the correlation time but the order parameter :math:`S^2_{int}` of the intermediate motion. 
    
.. _Razaei-Ghaleh et al. (2018): https://doi.org/10.1002/anie.201808172

The resulting spectral density function is then:

.. math::

    J(\omega) = S^2_{slow} \frac{\tau_{slow}}{1 + (\omega \tau_{slow})^2} + (1 - S^2_{slow}) S^2_{int} \frac{\tau_{int}}{1 + (\omega \tau_{int})^2}


The remaining **TODO the remaining what???** being accounted for by a fast motion with short (ps) correlation time (equation 18 in `Salvi et al. (2017)`_).

.. _Salvi et al. (2017): https://www.sciencedirect.com/science/article/pii/S0079656517300213

Given that

.. math::

    c = \dfrac{\eta_{xy}}{(3 cos^2 \theta - 1) p d} = 4 J(0) + 3 J(\omega_N)

where :math:`\eta_{xy}` is the cross-correlated relaxation rate, :math:`\theta` is the angle between the chemical shift anisotropy (CSA) tensor and the N-H bond vector, :math:`p` is the dipolar coupling constant, :math:`d` is the chemical shielding anisotropy constant, and :math:`\omega_N` is the Larmor frequency of the nitrogen nucleus. 
We can extract :math:`S^2_{int}` from the obtained value of :math:`c`, assuming that :math:`S^2_{slow}`, :math:`\tau_{slow}`, and :math:`\tau_{int}` are known:

.. math::

    S^2_{int} = \frac{|\eta_{xy}/(d*c*B_0*P_2(\cos(\theta)))| - 4(S^2_{slow} J(0, \tau_{slow}) + (1-S^2_{slow}) J(0, \tau_{fast})) - 3(S^2_{slow} J(\omega_N, \tau_{slow}) + (1-S^2_{slow}) J(\omega_N, \tau_{fast}))}{4(1-S^2_{slow})(J(0, \tau_{int}) - J(0, \tau_{fast})) + 3(1-S^2_{slow})(J(\omega_N, \tau_{int}) - J(\omega_N, \tau_{fast}))}


In this case the call is, for instance for synuclein at 600 MHz:
::

    t1t2ne tract --basedir <basedir> --tract <TRACT experiment number> --idp --MW 14.4 --selectregion --plot

The results are shown in the following figure.

.. figure:: _static/tract_synuclein.png
    :name: tract_synuclein
    :width: 80.0%

    The result of the TRACT analysis for the IDP synuclein at 600 MHz in the ``--selectregion`` mode. The raw data are provided in the ``examples`` directory.

At the end of the analysis, the software provides a command to generate the lists for running T\ :sub:`1`\  and T\ :sub:`2`\  experiments based on the obtained correlation time.

.. rubric:: Examples

**TODO** da popolare



.. rubric:: TODO

Non si capisce dove parli di IDP e dove no. Popola la sezione *Esempi* spostando i pezzi dove invochi il codice e spieghi a cosa serve cosa. Nel testo principale lascia la teoria e i risultati generali.


