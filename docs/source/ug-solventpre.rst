.. _sec;solventpre:

The ``solventpre`` subcommand
*****************************

This subcommand generates the lists for running :math:`^1` H T\ :sub:`1`\  and T\ :sub:`2`\  solvent PRE experiments. The expected relaxation rates are computed with the Freed Outer Sphere model. 
In this model, the spectral density function is computed as:

.. math::

    \tau_D = \frac{d^2}{D_{target} + D_{cosolute}} \\
    z = \sqrt{2 |\omega| \tau_D + \frac{\tau_D}{\tau_1}} \\
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

It can be called providing only the field as:

::

    t1t2ne solventpre --Larmor 600

Alternatively, the user can provide the expected intrinsic T\ :sub:`1`\  of the system and/or the expected linewidth:

::

    t1t2ne solventpre --Larmor 600 --T1 1 --lw 10

The expected linewidth can also be estimated from the molecular weight of the system:

::

    t1t2ne solventpre --Larmor 600 --MW 8.6

By default, the software will suggest a list of 8 delays for T\ :sub:`1`\  and 2 points for T\ :sub:`2`\ , as described by Iwahara, Tang, and Clore in `10.1016/j.jmr.2006.10.003`_. The number of points can be changed with the ``--nT`` argument.
    
    .. _10.1016/j.jmr.2006.10.003: https://doi.org/10.1016/j.jmr.2006.10.003


**Important**: The delays are calculated under the assumption that the pulse sequence used for T\ :sub:`2`\  is the one found in figure 1 of `10.1016/j.jmr.2006.10.003`_.:
