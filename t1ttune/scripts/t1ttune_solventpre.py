#! /usr/bin/env python3

import numpy as np

from .base import BaseCommand
from .textcolor import textcolor
from . import t1ttune_utils, f_ParaMeters_relax
import klassez as kz
import random

class SPREListsCmd(BaseCommand):
    SHORT_HELP = "Setup the time increment for solvent PRE experiments"
    DESCRIPTION = '\n'.join([
    'With this module, the time increment for solvent PRE experiments are computed.' 
    ,])
    
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--MW', type=float, help='The molecular weight of the protein in kDa, to be used for estimating the linewidth.')
        parser.add_argument('--lw', type=float, default=10, help='The expected linewidth in Hz, to be used for estimating the tau_slow in the IDP model. Default is 10 Hz.')
        parser.add_argument('--T1', type=float, default=1, help='The expected T1 in seconds, to be used for estimating the tau_slow in the IDP model. If not provided, it will be estimated from the molecular weight if provided, otherwise a default value of 1 s will be used.')
        parser.add_argument('--T', type=float, help='The Temperature in Kelvin')
        parser.add_argument('--nT', nargs='*', help='The number of increments in the suggested vdlist. For this module only one value must be provided. Default is 8.')
        parser.add_argument('--r', type=float, default=3.6, help='The distance of closest approach in Angstroms. Default is 3.6e-10 m, corresponding to the sum of the ionic radius of Gd^3+ and the van der Waals radius of water.')
        parser.add_argument('--B0', nargs='*', help='The magnetic field strength in Tesla. If not provided, the script will try to load it from the config file. If it is not found in the config file, an error will be raised.')
        parser.add_argument('--Larmor', nargs=1, type=float, help='The Larmor frequency in MHz, alternative to --B0.')
        parser.add_argument('--nucs', nargs='*', help='The nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.')
        parser.add_argument('--taue', type=float, help='The electron relaxation time in seconds. If not provided, a default value of 1e-9 s will be used.')
        parser.add_argument('--tauv', type=float, help='The correlation time for transient ZFS fluctuations in seconds. If not provided, a default value of 2.6e-11 s will be used, corresponding to Gd^3+ complexes.')
        parser.add_argument('--deltat', type=float, default=0.014, help='The transient ZFS parameter in cm^-1. Default is 0.014 cm^-1 for Gd^3+ complexes, None for no transient ZFS.')
        parser.add_argument('--c', type=float, help='The concentration of the paramagnetic cosolute in mM. If not provided, a default value of 1 mM will be used.')
        parser.add_argument('--access', type=float, help='The fraction of the sphere of accessibility for the cosolute. If not provided, ta default value of 0.5 will be used.')
        parser.add_argument('--D', nargs='*', type=float, help='A list of two diffusion coefficients in m^2/s, the first one for the target molecule and the second one for the cosolute. If not provided, default values of 1e-10 m^2/s for the target and 2.6e-10 m^2/s for the cosolute will be used.')
        parser.add_argument('--S', type=float, help='The spin quantum number of the paramagnetic center. If not provided, a default value of 3.5 will be used, corresponding to Gd^3+.')
        parser.add_argument('--AMe', type=float, help='The hyperfine coupling constant to the metal center in Hz. If not provided, a default value of None will be used, not used for Gd-DOTA.')
        parser.add_argument('--I', type=float, help='The nuclear spin quantum number. If not provided, a default value of None will be used, not used for Gd-DOTA.')
        parser.add_argument('--g', type=float, help='The electron g-factor. If not provided, a default value of None will be used.')
    @staticmethod
    def run(args):
        CO = t1ttune_utils.Conf_Optns(args, module='solventpre')
 
        solventpre(CO)
        t1ttune_utils.the_end(CO)    
        exit()        

def solventpre(CO):
    r"""
    Suggests a vdlist for the solventPRE experiments based on the Outer sphere model. The vdlist is geometrically spaced from 20 us to the longest delay. The user can choose the number of increments in the vdlist. The function prints the suggested vdlist in Bruker-readable format and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.

    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to load the data and perform the analysis.
        It is expected that the ``CO`` object has the following attributes:
        - ``B_0``: the magnetic field strength in Tesla
        - ``taue``: the electron relaxation time in seconds
        - ``tauv``: the correlation time for transient ZFS fluctuations in seconds default 2.6e-11 s, for Gd^3+ complexes)
        - ``c``: the concentration of the paramagnetic cosolute in mM. Default is 1 mM.
        - ``access``: the fraction of the sphere of accessibility for the cosolute. Default is 0.5.
        - ``D``: a list of two diffusion coefficients in m^2/s, the first one for the target molecule and the second one for the cosolute. If not provided, default values of 1e-10 m^2/s for the target and 2.6e-10 m^2/s for the cosolute are used.
        - ``S``: the spin quantum number of the paramagnetic center. Default is 3.5, corresponding to Gd^3+.
        - ``nucs``: the nuclei to use for the calculation of the relaxation rates. Default is ['1H', '15N'].  
        - ``r``: the distance of closest approach in meters. Default is 3.6e-10 m.
        - ``AMe``: the hyperfine coupling constant to the metal center in Hz. Default is None, not used for Gd-DOTA.
        - ``I``: the nuclear spin quantum number. Default is None, not used for Gd-DOTA.
        - ``g``: the electron g-factor. Default is None, uses free electron g-factor, not used for Gd-DOTA.
        - ``nT``: the number of increments in the suggested vdlist. Default is 8.
    
    Returns
    -------
    None
         The function prints the suggested vdlist in Bruker-readable format and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.

    """
    if hasattr(CO, 'nT'):
        nT1 = CO.nT[0] if isinstance(CO.nT, list) else CO.nT
        nT2 = CO.nT[1] if isinstance(CO.nT, list) and len(CO.nT) > 1 else CO.nT
    else:
        nT1 = 8
        nT2 = 2
    CO.add_ref('bertini')
    if not hasattr(CO, 'taue') and not hasattr(CO, 'tauv'):
        taue = 1e-9
        tauv = 2.6e-11
    if hasattr(CO, 'taue') and CO.taue is not None:
        taue = CO.taue
        tauv = None
    elif hasattr(CO, 'tauv') and CO.tauv is not None:
        tauv = CO.tauv
        taue = None
    else:
        taue = 1e-9
        tauv = 2.6e-11


    Gamma1, Gamma2 = f_ParaMeters_relax.OuterSphere(CO.B_0, c=CO.c, d=CO.r, D_target = CO.D[0], D_cosolute = CO.D[1], f=CO.access, taue=taue, tauv=tauv, deltat=CO.deltat, AMe=CO.AMe, I=CO.I, g=CO.g, S=CO.S, nuc=CO.nucs[0])


    if hasattr(CO, 'tau') and CO.tau is not None:
        tau_average = CO.tau[0]
        R2H = tau_average*7.42e9
        amidelinewidth = R2H/(np.pi)    
    else:
        if hasattr(CO, 'lw') and CO.lw is not None:
            amidelinewidth = CO.lw
        else:
            amidelinewidth = 10  # Default value
    T1 = CO.T1 if hasattr(CO, 'T1') and CO.T1 is not None else 1
    R1 = Gamma1 + T1
    R2 = Gamma2 + amidelinewidth/np.pi
    print(textcolor(f"\nEstimated R1 at the closest approach distance: {R1:.2f} s^-1", "blue"))
    print(textcolor(f"Estimated R2 at the closest approach distance: {R2:.2f} s^-1", "blue"))
    vdlist_PRE = np.geomspace(2e-5, 2/R1, num=nT1) #geometrically spaced list from 20us to 2*tau_average
    if CO.options['randomize']:
        random.shuffle(vdlist_PRE)
    print(textcolor('\nSuggested vdlist for T1 experiment:', 'blue'))
    print('-'*38)
    t1ttune_utils.out_vdlist(vdlist_PRE)
    if nT2 == 2:
        print(textcolor(f"Optimal Tb-Ta difference for optimizing Gamma2 measurement: {1.15/(R2):.5f} s", "blue"))
    else:
        vdlist_PRE_R2 = np.geomspace(2e-5, 2/(R2), num=nT2)
        if CO.options['randomize']:
            random.shuffle(vdlist_PRE_R2)
        print(textcolor('\nSuggested vdlist for R2 experiment:', 'blue'))
        print('-'*38)
        t1ttune_utils.out_vdlist(vdlist_PRE_R2) 
    
