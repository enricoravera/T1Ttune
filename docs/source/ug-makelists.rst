.. _sec;makelists:

The ``makelists`` subcommand
****************************

This subcommand generates the lists for running :math:`^15`N T\ :sub:`1`\  and T\ :sub:`2`\  relaxation experiments given the correlation times or the molecular weight.
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

