#! /usr/bin/env python3

"""
T1Ttune package for computing the optimal parameters for protein dynamics NMR experiments.

The package includes modules for computing:

    * Relaxation rates (R1, R2, NOE) based on a spectral density function and the provided correlation times and order parameters.
    * The optimal set of relaxation delays for each experiment (module `t1ttune_makelists`).
    * The optimal number of scans for each experiment (module `t1ttune_ns`).
    * The optimal list for setting up a TRACT experiment (module `t1ttune_setuptract`).
    * The correlation times or order parameters from a TRACT experiment (module `t1ttune_tract`).
            
"""

__version__ = "0.0.1a"
__author__ = "Enrico Ravera, Francesco Bruno, Letizia Fiorucci"
