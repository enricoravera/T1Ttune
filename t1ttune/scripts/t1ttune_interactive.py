#! /usr/bin/env python3

import os
import numpy as np

from .base import BaseCommand
from .textcolor import textcolor
import klassez as kz
from . import t1ttune_utils, fun_hetrelax_models, f_findfs
import random

class InteractiveCmd(BaseCommand):
    SHORT_HELP = "Interactive setup of the experiment"
    DESCRIPTION = '\n'.join([
    'With this module, the vdlist and vclist for protein dynamics experiment is computed.',
    'It can run without any input, in which case it will start from a default of 1 s for T1 and 3 loops for the CPMG block, and then suggest a vdlist based on the estimated R1_ave and R2_ave of the system.', 
    'In the standard operation mode, the tau_c or MW are used to compute the optimal vdlist.',
    'S2 is assumed to be 0.9 for the standard operation mode.',
    'The IDP mode requires input of the molecular weight of the protein to estimate the correlation time of the slow motion.',
    'Optionally, the length of the correlation window can be provided as an estimate for the intermediate correlation time.',
    'If not provided, it assumes a correlation time of 1.6 ns for the intermediate motion.'
    'Two S2 values are expected in the IDP mode. If not provided, they are assumed to be 0.15 for the slow motion and 0.35 for the intermediate motion.'
    ,])
    
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--basedir', type=str, help='The base directory for the experiment. Default is the current directory.')
        parser.add_argument('--xred', nargs='*', help='The xred values to use for the calculation of the lists.')
        parser.add_argument('--S2', nargs='*', help='The Lipari-Szabo order parameter S2 to use for the calculation of tau_c. In IDP mode, two values should be provide, else only one')
        parser.add_argument('--idp', action='store_true', help='Whether to use the IDP model to extract the order parameter S2 instead of tau_c.')
        parser.add_argument('--MW', type=float, help='The molecular weight of the protein in kDa, to be used for estimating the tau_slow in the IDP model. If not specified, will raise an error if --idp is used.')
        parser.add_argument('--corr_window_idp', type=int, default=20, help='To estimate the correlation time of the intermediate motion in the IDP model, a correlation time of a peptide of 20 residues is used by default. This parameter allows to change this value if needed.')
        parser.add_argument('--T', type=float, default=298.15, help='The temperature in Kelvin, to be used for estimating the tau_slow in the IDP model. Default is 298 K.')
        parser.add_argument('--nT', nargs='*', help='The number of increments in the suggested vdlist. For this module only one value must be provided. Default is 8.')
        parser.add_argument('--r', type=float, default=1.02, help='The length of the 1H-15N bond in Angstroms. Default is 1.02 A.')
        parser.add_argument('--Deltasigma', type=float, default=-160, help='The chemical shift anisotropy of the 15N nucleus in ppm. Default is -160 ppm.')
        parser.add_argument('--theta', type=float, default=17, help='The angle between the 1H-15N bond and the principal axis of the CSA tensor in degrees. Default is 17 degrees.')
        parser.add_argument('--B0', nargs='*', help='The magnetic field strength in Tesla. If not provided, the script will try to load it from the config file. If it is not found in the config file, an error will be raised.')
        parser.add_argument('--Larmor', type=float, help='The Larmor frequency of the nucleus in MHz. If not provided, the script will try to load it from the config file. If it is not found in the config file, an error will be raised.')
        parser.add_argument('--logscale', action='store_true', help='Whether to use a logarithmic scale for the vdlist and vclist. Default is False, which means a linear scale will be used.')
        parser.add_argument('--nucs', nargs='*', default=['1H', '15N'], help='The nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.')
        parser.add_argument('--large', action='store_true', help='Whether to create the lists for the "large" sequence, which is optimized for short T2 times. If True, the d21 value is set to 450 us and only 8 cycles per CPMG block are used instead of 16. Default is False.')
        parser.add_argument('--small', action='store_true', help='Whether to create the lists for the ".idp" sequence, which is optimized for long T2 times. If True, the d21 value is set to 600 us. Default is False.')
        parser.add_argument('--randomize', action='store_true', help='Whether to randomize the order of the values in the lists. Default is False.')

    @staticmethod
    def run(args):
        CO = t1ttune_utils.Conf_Optns(args, module='interactive')
 
        interactive_setup(CO)
        t1ttune_utils.the_end(CO)    
        exit()        
        
def interactive_setup(CO):
    """
    Interactive setup of the experiment. This script can only run on a spectrometer.
    
        Parameters
        ----------
        CO: Conf_Optns object
            The configuration options object containing the parameters for the experiment.
        
        Returns
        -------
        None
            The function prints the suggested vdlist and vclist.      
    """
    CO.add_ref('ferrage')
    CO.add_ref('Fiorucci')
    CO.add_ref('klassez')    

    CO.get_B0(config_p=t1ttune_utils.load_config())

    if hasattr(CO, 'tau') and CO.tau is not None:
        CO.add_ref('fushman')
        R1, R2, nOe = fun_hetrelax_models.R1R2nOe(CO.B_0, r=CO.r, nuc1=CO.nucs[0], nuc2=CO.nucs[1], Deltasigma=CO.Deltasigma, func=fun_hetrelax_models.LS_iso, f_args=(CO.S2, CO.tau))
        print(textcolor('Based on the provided parameters, the estimated R1 is {:.2f} s^-1 and R2 is {:.2f} s^-1.'.format(R1, R2), 'green'))
        T1_30 = -np.log(0.3)/R1
        T2_30 = -np.log(0.3)/R2
    else:
        T1_30 = 0.5
        T2_30 = 0.065

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
    logscale = CO.options['logscale']
    nT1 = CO.nT[0]
    if len(CO.nT) == 1:
        nT2 = nT1
    else:
        nT2 = CO.nT[1]
    
    if CO.basedir is None:
        print('\nFrom now on you need to run this script on the spectrometer, as it will need to access the acquisition parameters to suggest the lists for the experiment.')
        info = f_findfs.find_topspin()
        if not info['found'] or not info['spectrometer']:
            textcolor('TopSpin spectrometer installation not found. Please make sure you are running this script on the spectrometer and that TopSpin is properly installed.', 'red', bold=True)
            t1ttune_utils.the_end(CO) 
        CO.basedir = input('Please enter the base directory for the experiment (default is current directory): ') or '.'
    
        
    print('\nStarting the interactive setup of the T1 experiment...')
    print('Acquire the reference spectrum first with d7 = 20u.')
    refno = input('Please enter the reference spectrum number (default is 1): ') or '1'

    path_ref = os.path.join(CO.basedir, refno)
    if not os.path.exists(path_ref):
        raise RuntimeError('The specified reference spectrum number does not exist in the base directory. Please make sure to acquire the reference spectrum and enter the correct number.')
    try:
        S_ref = kz.Spectrum_1D(path_ref)
    except Exception as e:
        raise RuntimeError(f'Error occurred while loading the reference spectrum: {e}')
    S_ref.procs['wf']['mode'] = 'em'
    S_ref.procs['wf']['lb'] = 1/S_ref.acqus['AQ']
    print(f'Line broadening applied: {S_ref.procs["wf"]["lb"]:.2g} Hz')
    # Zerofill to twice the size
    S_ref.procs['zf'] = 2 * S_ref.fid.shape[-1]
    # Apply and do FT
    S_ref.process()
    # Remove digital filter
    S_ref.pknl()
    S_ref.adjph()
    S_ref.integrate(filename='refspec')
    limits = kz.misc.key_to_limits(list(S_ref.integrals.keys()))

    ta = S_ref.ngdict['acqus']['D'][7]
    whilecontrol = True
    iteration = 0
    while whilecontrol:
        print(textcolor(f'To achieve a reduction of the signal intensity of about {(1-CO.T1red)*100:.0f}%, try setting d7 to ' + t1ttune_utils.f4(T1_30), 'blue'))
        expno = input('Please enter the experiment number for the T1 experiment once it is done (default is 2): ') or '2'
        path_t1 = os.path.join(CO.basedir, expno)
        if not os.path.exists(path_t1):
            raise RuntimeError('The specified T1 experiment number does not exist in the base directory. Please make sure to acquire the T1 experiment and enter the correct number.')
        try:
            S_t1 = kz.Spectrum_1D(path_t1)
        except Exception as e:
            raise RuntimeError(f'Error occurred while loading the T1 spectrum: {e}')
        S_t1.procs['wf']['mode'] = 'em'
        S_t1.procs['wf']['lb'] = 1/S_t1.acqus['AQ']
        S_t1.procs['zf'] = 2 * S_t1.fid.shape[-1]
        S_t1.process()
        S_t1.pknl()
        S_t1.adjph()
        S_t1.integrate(filename='t1spec', limits=limits)
        tb = S_t1.ngdict['acqus']['D'][7]
        ratio = 0
        for key in S_t1.integrals.keys():
            ratio += S_t1.integrals[key] / S_ref.integrals[key] * (1/len(S_t1.integrals.keys()))
        R1_ave = -np.log(ratio) / (tb-ta)
        print(f'Iteration {iteration}: Estimated R1 from the reference and T1 spectra: {R1_ave:.2f} s^-1')
        T1_30_est = -np.log(0.3)/R1_ave
        if abs(T1_30_est - T1_30) / T1_30 > 0.1:
            T1_30 = T1_30_est
        else:
            whilecontrol = False
        iteration += 1
    print(textcolor(f'Convergence achieved for T1 estimation. Estimated $R1_{{ave}}$ is {R1_ave:.2f} s.', 'green'))
    print(textcolor(f'Set the recovery delay for the hetnOe experiment to at least {6/R1_ave:.3f}.\n', 'blue', bold=False))       
    T1max = -np.log(CO.T1red)/R1_ave
    if logscale:
        vdlist_T1 = [2e-5 - 1/R1_ave * np.log(1-(1-CO.T1red)*i/(nT1-1)) for i in range(nT1)] #logarithmically spaced list from 20u to to T1max
    else:
        vdlist_T1 = np.linspace(2e-5, T1max, num=nT1) #linearly spaced list from 20u to to T1max
    print(textcolor('\nvdlist for T1 experiment:', 'blue'))
    print('-'*25)
    if CO.options['randomize']:
        random.shuffle(vdlist_T1)
    t1ttune_utils.out_vdlist(vdlist_T1)


    print('\nStarting the interactive setup of the T2 experiment...')
    print('Acquire the reference spectrum first with l1 = 0.')
    if CO.options['idp'] and not CO.options['small']:
        CO.options['small'] = True
    if CO.options['small']:
        print(textcolor('Using ".idp" sequence, which is optimized for long T2 times. d21 = 600u', 'blue'))
        d21 = 600
    else:
        d21 = float(input('Enter the d21 value in microseconds (default 450): ').strip() or "450")
    p30 = float(input('Enter the p30 value in microseconds (default 80): ').strip() or "80")
    d31 = (p30*16+d21*32)
    if CO.options['large']:
        d31 = d31/2

    
    refno = input('Please enter the reference spectrum number (default is 1): ') or '1'

    path_ref = os.path.join(CO.basedir, refno)
    if not os.path.exists(path_ref):
        raise RuntimeError('The specified reference spectrum number does not exist in the base directory. Please make sure to acquire the reference spectrum and enter the correct number.')
    try:
        S_ref = kz.Spectrum_1D(path_ref)
    except Exception as e:
        raise RuntimeError(f'Error occurred while loading the reference spectrum: {e}')
    S_ref.procs['wf']['mode'] = 'em'
    S_ref.procs['wf']['lb'] = 1/S_ref.acqus['AQ']
    print(f'Line broadening applied: {S_ref.procs["wf"]["lb"]:.2g} Hz')
    # Zerofill to twice the size
    S_ref.procs['zf'] = 2 * S_ref.fid.shape[-1]
    # Apply and do FT
    S_ref.process()
    # Remove digital filter
    S_ref.pknl()
    S_ref.adjph()
    S_ref.integrate(filename='refspec')
    limits = kz.misc.key_to_limits(list(S_ref.integrals.keys()))

    ta = S_ref.ngdict['acqus']['L'][1]*d31*1e-6
    whilecontrol = True
    iteration = 0
    l1_30 = int(T2_30/(d31*1e-6))
    if l1_30 <= 0:
        l1_30 = 1

    while whilecontrol:
        print(textcolor(f'To achieve a reduction of the signal intensity of about {(1-CO.T2red)*100:.0f}%, try setting l1 to {l1_30:.0f}'))
        expno = input('Please enter the experiment number for the T2 experiment once it is done (default is 2): ') or '2'
        path_t2 = os.path.join(CO.basedir, expno)
        if not os.path.exists(path_t2):
            raise RuntimeError('The specified T2 experiment number does not exist in the base directory. Please make sure to acquire the T2 experiment and enter the correct number.')
        try:
            S_t2 = kz.Spectrum_1D(path_t2)
        except Exception as e:
            raise RuntimeError(f'Error occurred while loading the T2 spectrum: {e}')
        S_t2.procs['wf']['mode'] = 'em'
        S_t2.procs['wf']['lb'] = 1/S_t2.acqus['AQ']
        S_t2.procs['zf'] = 2 * S_t2.fid.shape[-1]
        S_t2.process()
        S_t2.pknl()
        S_t2.adjph()
        S_t2.integrate(filename='t2spec', limits=limits)
        tb = S_t2.ngdict['acqus']['L'][1]*d31*1e-6
        ratio = 0
        for key in S_t2.integrals.keys():
            ratio += S_t2.integrals[key] / S_ref.integrals[key] * (1/len(S_t2.integrals.keys()))
        R2_ave = -np.log(ratio) / (tb-ta)
        print(f'Iteration {iteration}: Estimated R2 from the reference and T2 spectra: {R2_ave:.2f} s^-1')
        l1_30_est = int((-np.log(0.3)/R2_ave)/(d31*1e-6))
        if l1_30_est == 0:
            l1_30_est = 1
        if l1_30_est != l1_30:
            l1_30 = l1_30_est
        else:
            whilecontrol = False
        iteration += 1
    
    print(textcolor(f'Convergence achieved for T2 estimation. Estimated $R2_{{ave}}$ is {R2_ave:.2f} s.', 'green'))
    T2max = -np.log(CO.T2red)/R2_ave
    nmax = int(T2max/(d31*1e-6))
    if T2max > 0.250:
        print(textcolor(f'Alert: the longest CPMG block for a residual signal of {CO.T2red*100:.0f}% is {T2max:.2f} s, which is likely too long for any equipment.', 'red', bold=True))
        print('Setting the maximum CPMG block to to 250ms')
        T2max = 0.250
        if CO.options['large']:
            print(textcolor('You chose the "large" sequence, which is optimized for short T2 times, this option will be disabled.', 'yellow', bold=True))
            CO.options['large'] = False

    #print(f'd21 = {d21} us, p30 = {p30} us, cpmgblock = {d31} us')
    print('\n')
    print(textcolor(f'The longest CPMG block for T2 for a residual signal of {CO.T2red*100:.0f}% should be {T2max:.2f} s, with {nmax} loops.' , 'default', bold=True)) #grassetto
    print(textcolor('Check if this is too long for your equipment before running the experiment', 'red')) #rosso
    if nT2 > nmax:
        print(textcolor(f'Warning: the number of increments you chose for the CPMG experiment is {nT2}, which is more than the recommended number of loops {nmax} for a longest CPMG block of {T2max:.2f} s. The list will contain duplicates.', 'yellow')) #giallo
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

    