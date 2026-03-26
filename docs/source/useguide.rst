.. |nbsp| unicode:: U+00A0
   :trim:


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

