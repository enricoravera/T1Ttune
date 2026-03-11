#! /usr/bin/env python3

import numpy as np
import random

from .base import BaseCommand
from .textcolor import textcolor
from . import t1ttune_utils, fun_hetrelax_models

class MakeListsCmd(BaseCommand):
    """Make lists of values for the protein dynamics experiments"""
    SHORT_HELP = "Make lists of values for the protein dynamics experiments"
    DESCRIPTION = '\n'.join([
    'With this module, the lists of values for the protein dynamics experiments are created, and the heteronuclear NOE expected value is printed.',
    'In the standard operation mode, the tau_c or MW are used to compute the optimal lists.'
    'S2 is assumed to be 0.9 for the standard operation mode.'
    'The IDP mode requires input of the molecular weight of the protein to estimate the correlation time of the slow motion.',
    'Optionally, the length of the correlation window can be provided as an estimate for the intermediate correlation time.',
    'If not provided, it assumes a correlation time of 1.6 ns for the intermediate motion.'
    'Two S2 values are expected in the IDP mode. If not provided, they are assumed to be 0.15 for the slow motion and 0.35 for the intermediate motion.'
    ,])
    
    @staticmethod
    def add_arguments(parser):
        """Adds the arguments for the MakeLists command to the parser."""
        parser.add_argument('--S2', nargs='*', help='The Lipari-Szabo order parameter S2 to use for the calculation of tau_c. In IDP mode, two values should be provided, else only one')
        parser.add_argument('--tau', nargs='*', help='The correlation time tau to use for the calculation. In IDP mode, two values should be provided, else only one.')
        parser.add_argument('--idp', action='store_true', help='Whether to use the IDP model to extract the order parameter S2 instead of tau_c.')
        parser.add_argument('--MW', type=float, help='The molecular weight of the protein in kDa, to be used for estimating the tau_slow in the IDP model. If not specified, will raise an error if --idp is used.')
        parser.add_argument('--corr_window_idp', type=int, help='To estimate the correlation time of the intermediate motion in the IDP model, a correlation time of a peptide of 20 residues is used by default. This parameter allows to change this value if needed.')
        parser.add_argument('--T', type=float, default=298.15, help='The temperature in Kelvin, to be used for estimating the tau_slow in the IDP model. Default is 298 K.')
        parser.add_argument('--xred', nargs='*', help='The percent residual signal for the longest delay of the experiments. Accepts a list.')
        parser.add_argument('--nT', nargs='*', type=int, default=[8, 8], help='The number of increments in the suggested vdlist and vclist. Default is 8.')
        parser.add_argument('--r', type=float, default=1.02, help='The length of the 1H-15N bond in Angstroms. Default is 1.02 A.')
        parser.add_argument('--Deltasigma', type=float, default=-160, help='The chemical shift anisotropy of the 15N nucleus in ppm. Default is -160 ppm.')
        parser.add_argument('--theta', type=float, default=17, help='The angle between the 1H-15N bond and the principal axis of the CSA tensor in degrees. Default is 17 degrees.')
        parser.add_argument('--B0', nargs='*', type=float, help='The magnetic field strength in Tesla. If not provided, the script will try to load it from the config file. If it is not found in the config file, an error will be raised.')
        parser.add_argument('--Larmor', nargs=1, type=float, help='The Larmor frequency in MHz, alternative to --B0.')
        parser.add_argument('--nucs', nargs='*', default=['1H', '15N'], help='The nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.')
        parser.add_argument('--logscale', action='store_true', help='Whether to create the lists in logarithmic scale. If True, the protocol suggested by F. Ferrage is used. Default is False.')
        parser.add_argument('--large', action='store_true', help='Whether to create the lists for the "large" sequence, which is optimized for short T2 times. If True, the d21 value is set to 450 us and only 8 cycles per CPMG block are used instead of 16. Default is False.')
        parser.add_argument('--small', action='store_true', help='Whether to create the lists for the ".idp" sequence, which is optimized for long T2 times. If True, the d21 value is set to 600 us. Default is False.')
        parser.add_argument('--randomize', action='store_true', help='Whether to randomize the order of the values in the lists. Default is False.')
    @staticmethod
    def run(args):
        """Runs the MakeLists command."""
        CO = t1ttune_utils.Conf_Optns(args)
        R1, R2, nOe = fun_hetrelax_models.R1R2nOe(CO.B_0, r=CO.r, nuc1=CO.nucs[0], nuc2=CO.nucs[1], Deltasigma=CO.Deltasigma, func=fun_hetrelax_models.LS_iso, f_args=(CO.S2, CO.tau))
        CO.add_ref('fushman')
        print('\n')
        print(textcolor('Estimated values for heteronuclear relaxation:\n', 'green', bold=False))
        print(f'expected T1 = {1/R1:.3f} s\nexpected T2 = {1/R2:.3f} s\nexpected hetnOe = {nOe:.3f}')
        print(textcolor(f'Set the recovery delay for the hetnOe experiment to at least {6/R1:.3f}.\n', 'blue', bold=False))       
        create_lists(CO, R1, R2)
        t1ttune_utils.the_end(CO)
        exit()
        
def create_lists(CO, R1, R2):
    """
    Creates vdlist and vclist for T1 and T2 experiments based on the computed R1 and R2 values.
    The vdlist is created so that the longest delay of the T1 experiment corresponds to a residual signal intensity of T1red, which is a user-defined percentage reduction. Defaults to 5%.
    The vclist is created so that the longest CPMG block of the T2 experiment corresponds to a residual signal intensity of T2red, which is a user-defined percentage reduction. Defaults to 10%.
    The user can choose to create the lists in logarithmic or linear scale. In logscale, the protocol suggested by F. Ferrage in `his protocol`_ is used.
    
    .. _his protocol: https://doi.org/10.1007/978-1-61779-480-3_9
    
    
    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to load the data and perform the analysis. The following attributes of the CO object are used in this function:

        *   B_0 : float
            The magnetic field strength in Tesla.
        *   S2 : list of float
            The order parameters to use for the calculation. In IDP mode, two values should be provide, else only one.
        *   tau : list of float
            The correlation times to use for the calculations. In IDP mode, two values should be provide, else only one.
        *   nT : list of int
            The number of increments for the T1 and T2 experiments. nT[0] corresponds to T1 and nT[1] corresponds to T2. If only one value is provided, it is used for both
        *   xred: list of float
            The percent residual signal for the longest delay of the experiments. If only one value is provided, it is used for both T1 and T2. If not provided, defaults to 5% for T1 and 10% for T2.

        The following options of the CO object are also used in this function:

        *   large : bool, optional
            Whether to create the lists for the "large" sequence, which is optimized for short T2 times. If True, the d21 value is set to 450 us and only 8 cycles per CPMG block are used instead of 16. Default is False.
        *   small : bool, optional
            Whether to create the lists for the ".idp" sequence, which is optimized for long T2 times. If True, the d21 value is set to 600 us. Default is False.
        *   logscale : bool, optional
            Whether to create the lists in logarithmic scale. If True, the protocol suggested by F. Ferrage is used. Default is True.
        *   randomize : bool, optional
            Whether to randomize the order of the values in the lists. Default is False.
    R1 : float
        The R1 relaxation rate in s^-1.
    R2 : float
        The R2 relaxation rate in s^-1.
        
    Returns
    -------
    None
        The function prints the created vdlist and vclist in Bruker-readable format and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.
    
    """
    
    print(f'Creating vdlist with {CO.nT[0]} points and vclist with {CO.nT[1]} points...')
    
    logscale = CO.options['logscale']
    if logscale:
        CO.add_ref('ferrage')
    if hasattr(CO, 'xred'):
        if len(CO.xred) == 1:
            CO.T1red = CO.T2red = CO.xred[0]*0.01
        elif len(CO.xred) == 2:
            CO.T1red = CO.xred[0]*0.01
            CO.T2red = CO.xred[1]*0.01
        else:
            print(textcolor('Warning: the xred parameter accepts only one or two values, defaulting to 5% for T1 and 10% for T2', 'yellow', bold=True))        
    else:
        CO.T1red = float(input('Enter percent residual signal for the longest delay of the T1 experiment (default 5%): ').strip() or "5")*0.01
        CO.T2red = float(input('Enter percent residual signal for the longest CPMG block of the T2 experiment (default 10%): ').strip() or "10")*0.01        
    T1max = -np.log(CO.T1red)/R1
    T2max = -np.log(CO.T2red)/R2
    nT1 = CO.nT[0]
    if len(CO.nT) == 1:
        nT2 = nT1
    else:
        nT2 = CO.nT[1]
    if T2max > 0.250:
        print(textcolor(f'Alert: the longest CPMG block for a residual signal of {CO.T2red*100:.0f}% is {T2max:.2f} s, which is likely too long for any equipment.', 'red', bold=True))
        print('Setting the maximum CPMG block to to 250ms')
        T2max = 0.250
        if CO.options['large']:
            print(textcolor('You chose the "large" sequence, which is optimized for short T2 times, this option will be disabled.', 'yellow', bold=True))
            CO.options['large'] = False
    if CO.options['small']:
        print(textcolor('Using ".idp" sequence, which is optimized for long T2 times. d21 = 600u', 'blue'))
        d21 = 600
    else:
        d21 = float(input('Enter the d21 value in microseconds (default 450): ').strip() or "450")
    p30 = float(input('Enter the p30 value in microseconds (default 80): ').strip() or "80")
    d31 = (p30*16+d21*32)
    if CO.options['large']:
        d31 = d31/2
    else:
        if -np.log(1-(1-CO.T2red)*(1/(nT2-1))) < d31*1e-6*R2:
            print(textcolor('Warning', 'yellow') + f': the second CPMG point in logscale would be a repetition of the first, switching to "large" sequence') #giallo
            d31 = d31/2
    if CO.options['small']:
        if d31>0.5/R1:
            print(textcolor('Warning', 'yellow') + f': the CPMG block of {d31*1e-6:.2f} s is too long for the T2 timescale, switching off the  "small" option') #giallo
            d21 = 450
            d31 = (p30*16+d21*32)
            CO.options['small'] = False
        else:
            CO.add_ref('bolognesi')
        
    nmax = int(T2max/(d31*1e-6))
    #print(f'd21 = {d21} us, p30 = {p30} us, cpmgblock = {d31} us')
    print('\n')
    print(textcolor(f'The longest delay for the T1 experiment for a residual signal of {CO.T1red*100:.0f}% should be {T1max:.2f} s.', 'default', bold=True)) #grassetto
    print(textcolor(f'The longest CPMG block for T2 for a residual signal of {CO.T2red*100:.0f}% should be {T2max:.2f} s, with {nmax} loops.' , 'default', bold=True)) #grassetto
    print(textcolor('Check if this is too long for your equipment before running the experiment', 'red')) #rosso
    if nT2 > nmax:
        print(textcolor(f'Warning: the number of increments you chose for the CPMG experiment is {nT2}, which is more than the recommended number of loops {nmax} for a longest CPMG block of {T2max:.2f} s. The list will contain duplicates.', 'yellow')) #giallo
    #create the vdlist for the 15N T1 experiment
    if logscale:
        vdlist_T1 = [2e-5 - 1/R1 * np.log(1-(1-CO.T1red)*i/(nT1-1)) for i in range(nT1)] #logarithmically spaced list from 20u to to T1max
    else:
        vdlist_T1 = np.linspace(2e-5, T1max, num=nT1) #linearly spaced list from 20u to to T1max
    print(textcolor('\nvdlist for T1 experiment:', 'blue'))
    print('-'*25)
    if CO.options['randomize']:
        random.shuffle(vdlist_T1)
    t1ttune_utils.out_vdlist(vdlist_T1)
    #create the vclist for the 15N T2 experiment
    if logscale:
        vclist_T2 = [-1/R2 * np.log(1-(1-CO.T2red)*i/(nT2-1))/(d31*1e-6) for i in range(nT2)]
        vclist_T2 = np.array(vclist_T2, dtype=np.uint64) #convert to np.uint64
    else:
        vclist_T2 = np.linspace(0,nmax, num=int(nT2), dtype=np.uint64) #linearly spaced list from 0 to nmax
    vclist_T2 = list(vclist_T2)
    if CO.options['randomize']:
        random.shuffle(vclist_T2)
    print(textcolor('\nvclist for T2 experiment:', 'blue'))
    print('-'*25)
    for value in vclist_T2:
        print(f'{value:d}')
    print('\n')
