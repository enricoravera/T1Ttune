.. |nbsp| unicode:: U+00A0
   :trim:


##########
USER GUIDE
##########

************
Introduction
************

T1T2ne is an helper software package to setup protein dynamics :sup:`15`\ N T\ :sub:`1`\  and T\ :sub:`2`\ and hetNOE experiments. The idea is to save spectrometer time and to optimize the quality of the data, by automatizing the selection of the parametes.
The software is executed from command line, and the main entry point is the ``t1t2ne`` command. 
To get a list of available subcommands, you can run:
::

    t1t2ne --help

To get more information on a specific subcommand, you can run:
::
    
    t1t2ne <subcommand> --help

The following subcommands are available:


- Setup of relaxation experiments
    -   | ``interactive``: an interactive command line tools that guides the user through the selection of the appropriate parameters for :sup:`15`\ N relaxation measurements.
        | See section :ref:`sec;interactive`.
    -   | ``makelists``: a command line tool to generate the lists for running :sup:`15`\ N T\ :sub:`1`\  and T\ :sub:`2`\  relaxation experiments given the correlation times or the molecular weight.
        | See section :ref:`sec;makelists`.
    -   | ``ns``: a command line tool to estimate the number of scans required for a given experiment, starting from an HSQC experiment.
        | See section :ref:`sec;ns`.
- Setup and analysis of the TRACT experiment
    -   | ``setuptract``: a command line tool to set up a TRACT experiment on the spectrometer, providing the appropriate delay list based on the correlation time or molecular weight.
        | See section :ref:`sec;setuptract`.
    -   | ``tract``: a command line tool to fit a TRACT experiment and extract the average correlation time of the system.
        | See section :ref:`sec;tract`.
- Solvent PRE
    -   | ``solventpre``: a command line tool to generate the lists for running :sup:`1`\ H T\ :sub:`1`\  and T\ :sub:`2`\  solvent PRE experiments.
        | See section :ref:`sec;solventpre`.


Minimal Theoretical Background
******************************

We here provide the minimal theoretical background to understand the rationale behind the software. For a more detailed introduction we will provide in-line links to appropriate references.
Nuclear spin relaxation processes are caused by the fluctuation of the interactions of the nuclear spin with its environment. Such fluctuations can occur on a range of timescales, hence contribute to different frequencies.
This is described through the spectral density functions :math:`J(\omega)`, which are the Fourier transform of the correlation function of the fluctuating interactions.
For proteins, it is common to assume that the spectral density function can be described through the `Lipari-Szabo`_ model-free approach, which assumes that the motion of the N-H bond vector can be described through a global reorientation (slow motion) and an internal motion (fast motion).
The global reorientation is described through a correlation time :math:`\tau_c` and an order parameter :math:`S^2`, which describes the amplitude of the motion. The internal motion is described through an order parameter :math:`S^2_f = 1 - S^2` and a correlation time :math:`\tau_f`. The resulting spectral density function is then:

.. math::

    J(\omega) = \frac{2}{5} \left[ S^2 \frac{\tau_c}{1 + (\omega \tau_c)^2} + (1 - S^2) \frac{\tau_f}{1 + (\omega \tau_f)^2} \right]


which, assuming that :math:`\tau_f` is of the order of picoseconds, and thus :math:`\tau_f \ll \tau_c`, can be approximated to:

.. math::

    J(\omega) \approx \frac{2}{5} S^2 \frac{\tau_c}{1 + (\omega \tau_c)^2} 

For Intrinsically Disordered Proteins, several different motional models have been proposed (see `Salvi`_ for a review). In this software we have opted for the model proposed by `Parigi et al.`_, which assumes that the motion can be described through a slow motion, an intermediate motion and a fast motion. The slow motion is assumed to correspond to the tumbling of a protein of the same size of the one under study, and the order parameter of the slow motion is assumed to be small (e.g., 0.15). The intermediate motion is assumed to have a correlation time of a peptide of 20 residues (1.6 ns). The resulting spectral density function is then:

.. math::

    J(\omega) = \frac{2}{5} \left[ S^2_{slow} \frac{\tau_{slow}}{1 + (\omega \tau_{slow})^2} + (1 - S^2_{slow}) S^2_{int} \frac{\tau_{int}}{1 + (\omega \tau_{int})^2} + (1 - S^2_{slow}) (1 - S^2_{int}) \frac{\tau_{fast}}{1 + (\omega \tau_{fast})^2} \right]

For a :sup:`15`\ N-:sup:`1`\ H pair, the relevant interactions are the dipolar coupling and the chemical shielding anisotropy.
The dipolar coupling is described through the constants :math:`p` and :math:`\delta`.

.. math:: p = \frac{\mu_0 h \gamma_H \gamma_N}{16 \pi^2 r^3} 

.. math:: \delta = -\frac{\Delta \sigma \gamma_N}{3} 
    
 Where:

    -   :math:`\mu_0` is the vacuum permeability in T m/A
    -   :math:`h` is the Planck constant in J s
    -   :math:`\gamma_H` is the gyromagnetic ratio of the proton in rad/s/T
    -   :math:`\gamma_N` is the gyromagnetic ratio of the nitrogen in rad/s/T
    -   :math:`r` is the distance between the nuclei in meters
    -   :math:`\Delta \sigma` is the anisotropy of the chemical shielding tensor, which is a dimensionless quantity. For a N-H pair, it is around 160 ppm.
    
The relaxation rates :math:`R_1`, :math:`R_2` and NOE enhancement factor :math:`\eta` at a given the magnetic field strength :math:`B` can be calculated as reported by `Fushman`_:


.. math:: R_1 = p^2 (J(\omega_H - \omega_N) + 6 J(\omega_H + \omega_N)) + 3 (\delta^2 B^2 + p^2) J(\omega_N)

.. math:: R_2 = R_{ex} + \frac{1}{2} p^2 (J(\omega_H - \omega_N) + 6 J(\omega_H + \omega_N) + 6 J(\omega_H)) + \frac{1}{2} (\delta^2 B^2 + p^2) (4 J(0) + 3 J(\omega_N))

.. math:: \eta = 1 - \frac{p^2 \gamma_H}{\gamma_N} \frac{6 J(\omega_H + \omega_N) - J(\omega_H - \omega_N)}{R_1}

.. note::

    Throughout the code, it is assumed that :math:`R_{ex}` is negligible and that also the contribution from neighboring protons is currently neglected.

Given the computed values of :math:`R_1`, :math:`R_2` and :math:`\eta`, it is possible to compute the an optimal sampling scheme for the experiments, as described by `Ferrage`_.

The use of the symbol :math:`\eta` is quite unfortunate, as it also indicates the cross-correlated relaxation rates.
For the cross-correlation between the dipolar coupling and the CSA, the rates are reported by `Salvi`_:
   
.. math:: \eta_z = \frac{1}{15} P_2(\cos\theta) p \delta B 6 J(\omega_X)

.. math:: \eta_{xy} = \frac{1}{15} P_2(\cos\theta) p \delta B (4 J(0) + 3 J(\omega_X))

where :math:`\theta` is the angle between the dipolar coupling and the CSA tensors, which for a N-H pair is around 17 degrees.

`Goldman`_ noted that:

.. math:: \eta_{red} = \dfrac{\eta_{xy}}{(3 cos^2 \theta - 1) p d} = 4 J(0) + 3 J(\omega_N)

From this relation and from the form of the spectral density function, it is possible to extract the correlation time of the system from the value of :math:`\eta_{xy}`. This is the rationale behind the TRACT experiment (`Lee et al. (2006)`_).
Recently, `Robson et al. (2021)`_ have provided a closed-form expression to extract the correlation time from the value of :math:`\eta_{xy}` including the possibility to specify a value for the order parameter. The equation is reported in the original paper and is too long to be reproduced here. 
We note that, while usually not advised, the same reasoning can be applied to extract the order parameter of an IDP:

.. math::

    S^2_{int} = \frac{\eta_{red} - 4(S^2_{slow} J(0, \tau_{slow}) + (1-S^2_{slow}) J(0, \tau_{fast})) - 3(S^2_{slow} J(\omega_N, \tau_{slow}) + (1-S^2_{slow}) J(\omega_N, \tau_{fast}))}{4(1-S^2_{slow})(J(0, \tau_{int}) - J(0, \tau_{fast})) + 3(1-S^2_{slow})(J(\omega_N, \tau_{int}) - J(\omega_N, \tau_{fast}))}

This closes the circle and allows to have a more robust estimate of the relaxation times and hence to obtain more meaningful delay lists, which in turn gives better experimental data. 

Solvent PRE and the Freed Outer Sphere Model
********************************************

`Solvent PRE`_ experiments are based on the fact that the relaxation of a nucleus can be enhanced by the presence of a paramagnetic agent in solution. The effect is distance-dependent, and can be used to obtain information on the solvent accessibility of different residues in a protein.
he expected relaxation rates are computed with the Freed Outer Sphere model. 
In this model, the spectral density function is computed as:

.. math::

    \tau_D &= \frac{d^2}{D_{target} + D_{cosolute}} \\
    z &= \sqrt{2 |\omega| \tau_D + \frac{\tau_D}{\tau_1}} \\
    J_{outer}(\omega) &= \frac{2}{5} \frac{1 + 5 z / 8 + z^2 / 8}{1 + z + z^2/2 + z^3/6 + 4 z^4/81 + z^5/81 + z^6/648}

.. FOR reference J_{outer}(\omega) &= \frac{2}{5} \frac{1 + {5 z}{8} + \frac{z^2}{8}}{1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{4 z^4}{81} + \frac{z^5}{81} + \frac{z^6}{648}}
        
where: 
    -   :math:`d` is the distance of closest approach of the paramagnetic center and the nucleus in meters
    -   :math:`D_{target}` is the diffusion coefficient of the target molecule in m^2/s
    -   :math:`D_{cosolute}` is the diffusion coefficient of the cosolute molecule in m^2/s
    -   :math:`\tau_1` is the correlation time for electron relaxation in seconds
    -   :math:`\omega` is the angular frequency where the spectral density is evaluated, in rad/s

Equations 6.42, 6.48 and 6.50 in `Bertini et al. 2016`_.
 

This spectral density function is then used to compute the expected PRE rates, with optional transient zero-field splitting (ZFS) contribution to the electron relaxation.
Default values are for 1 mM Gd-DOTA and a protein of 10 kDa at room temperature. The default values correspond to Gd-DOTA, and are taken from `Li et al. 2002`_.

.. math:: 

    k_{outer} = \frac{16\pi}{81} \left( \frac{\mu_0}{4\pi} \right)^2 h^2 \gamma_n^2 \gamma_e^2 S(S+1) f N_A \frac{c}{d(D_{target}+D_{cosolute})}\\
    R_1 = k_{outer} [7 J_{outer} (\omega_S) + 3 J_{outer} (\omega_I)]\\
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

    R_{1e} = \frac{2}{15} \Delta_t^2 \left(4S(S+1) - 3\right) \left(J(\omega_e) + 2J(2\omega_e)\right)\\
    R_{2e} = \frac{1}{15} \Delta_t^2 \left(4S(S+1) - 3\right) \left(3J(0) + 5J(\omega_e) + J(2\omega_e)\right)

Where:
    -   :math:`\Delta_t` is the transient zero-field splitting, provided in inverse centimeters and translated automatically to Hz
    -   :math:`\omega_e` is the electron Larmor frequency in rad/s


    .. _Lipari-Szabo: https://pubs.acs.org/doi/10.1021/ja00381a009
    .. _Parigi et al.: https://doi.org/10.1021/ja506820r
    .. _Ferrage: https://doi.org/10.1007/978-1-61779-480-3_9 
    .. _Fushman: https://pmc.ncbi.nlm.nih.gov/articles/PMC4361738/
    .. _Salvi: https://www.sciencedirect.com/science/article/pii/S0079656517300213
    .. _Goldman: https://doi.org/10.1016/0022-2364(84)90055-6
    .. _Lee et al. (2006): https://doi.org/10.1016/j.jmr.2005.08.014
    .. _Robson et al. (2021): https://doi.org/10.1007/s10858-021-00379-5
    .. _Solvent PRE: https://doi.org/10.1016/j.pnmrs.2022.09.001
    .. _Bertini et al. 2016: https://www.sciencedirect.com/science/chapter/monograph/pii/B9780444634368000065
    .. _Li et al. 2002: https://pubs.acs.org/doi/full/10.1021/ic0200390

*******************************
Setup of relaxation experiments
*******************************

.. include:: ug-interactive.rst
.. include:: ug-makelists.rst
.. include:: ug-ns.rst


******************************************
Setup and analysis of the TRACT experiment
******************************************

.. include:: ug-setuptract.rst
.. include:: ug-tract.rst


***********
Solvent PRE
***********

.. include:: ug-solventpre.rst


.. include:: acknowledgements.rst

