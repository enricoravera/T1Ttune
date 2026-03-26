.. _sec;setuptract:

The ``setuptract`` subcommand
*****************************

This subcommand generates the delay list for running a TRACT experiment on the spectrometer, providing the appropriate delay list based on the correlation time or molecular weight.
The general syntax for this command is the following:
::

    t1t2ne setuptract --MW <MW> --Larmor <Frequency in MHz>

The user needs to provide either the molecular weight or the correlation time of the system, and the Larmor frequency of \ :sup:`1`\ H at the given field. 
If the program is run on the spectrometer, this latter information is not needed, as it is direcly taken from there.

The software will then generate a delay list for the TRACT experiment that covers a range of delays that is appropriate for the expected correlation time, and that maximizes the `Fisher Information criterium`_ for matching the decay rates of both :math:`\alpha` and :math:`\beta` components of the TROSY.

.. seealso::

   :class:`t1t2ne.scripts.t1t2ne_setuptract`


.. rubric:: Examples

If I have a 8.6 kDa protein (ubiquitin) and a 600 MHz spectrometer, the script should be invoked as:

-   *Case 1:* I know the Molecular Weight and I am at the spectrometer computer
    ::

        t1t2ne setuptract --MW 8.6 

-   *Case 2:* I know the correlation time and I am at the spectrometer computer
    ::

        t1t2ne setuptract --tauc 666

-   *Case 3:* I know the Molecular Weight and I am `not` at the spectrometer computer
    ::

        t1t2ne setuptract --MW 8.6 --Larmor 600

-   *Case 4:* I know the correlation time and I am `not` at the spectrometer computer
    ::

        t1t2ne setuptract --tauc 666 --Larmor 600


.. rubric:: **TODO**

Qui bisogna popolare la sezione "esempi" coprendo i vari casi d'uso --> so il peso molecolare, so il tau, sono su uno spettrometro o no


