#! /usr/bin/env python3

import numpy as np
import argparse
import klassez as kz
import os

from . import t1t2ne_utils, f_fit, hydrodynamics_utils

from .base import BaseCommand
from .textcolor import textcolor
from .tract_extra import split_tract


class NSCmd(BaseCommand):
    SHORT_HELP = "Helps you select a reasonable number of scans for the dynamics experiments."
    DESCRIPTION = '\n'.join([
    'With this module, you can estimate the optimal number of scans for your NMR experiments based on the desired signal-to-noise ratio (SNR) and relaxation properties of the sample.',
    'The module reads a reference spectrum, then computes its SNR, and finally estimates the number of scans needed to achieve a target SNR for the T1, T2, and heteronuclear NOE experiments.',
    ])
    
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--basedir', type=str, help='Base directory for the experiment (where the reference spectrum is located)')
        parser.add_argument('--tract', type=str, help='expno of the TRACT experiment (used to estimate the SNR). Either this or --hsqc should be provided')
        parser.add_argument('--hsqc', type=str, help='expno of the reference HSQC spectrum (used to estimate the SNR). Either this or --tract should be provided')
        parser.add_argument('--nres', type=int, help='The number of non-proline residues in the protein. If not provided, it will be estimated from the molecular weight.')
        parser.add_argument('--lw', type=float, help='The linewidth of the peaks in the reference spectrum in Hz. If not provided, it will be estimated from the reference spectrum.')
        parser.add_argument('--xred', nargs='*', help='The percent residual signal for the longest delay of the experiment, or NOE efficiency as percentage. Accepts a list, computes the optimal number of scans for each value')
        parser.add_argument('--T', nargs=1, help='Temperature in Kelvin (used to recompute tau_c based on temperature)')
        parser.add_argument('--phase', '-p', action='store_true', help='Whether to phase the spectra.')
        parser.add_argument('--S2', nargs='*', help='The order parameters to use for the calculation of tau_c. In IDP mode, two values should be provide, else only one')
        parser.add_argument('--idp', action='store_true', help='Whether to use the IDP model to compute the spectral density function.')
        parser.add_argument('--MW', type=float, help='The molecular weight of the protein in kDa, used to estimate the number of residues of the protein. In the IDP model it is also used to compute tau_slow. If not specified, will raise an error if --idp is used.')
        parser.add_argument('--corr_window_idp', type=float, default=20, help='To estimate the correlation time of the intermediate motion in the IDP model, a correlation time of a peptide of 20 residues is used by default. This parameter allows to change this value if needed.')
        parser.add_argument('--r', type=float, default=1.02, help='The length of the 1H-15N bond in Angstroms. Default is 1.02 A.')
        parser.add_argument('--Deltasigma', type=float, default=-160, help='The chemical shift anisotropy of the 15N nucleus in ppm. Default is -160 ppm.')
        parser.add_argument('--theta', type=float, default=17, help='The angle between the 1H-15N bond and the principal axis of the CSA tensor in degrees. Default is 17 degrees.')
    @staticmethod
    def run(args: argparse.Namespace) -> None:
        CO = t1t2ne_utils.Conf_Optns(args, module='NS')
        SNR = estimate_snr(CO)
        suggest_scans(CO, SNR)
        t1t2ne_utils.the_end(CO)
        
        
def estimate_snr(CO):
    """
    Provides an estimate of the expected signal-to-noise ratio (SNR) per scan for a single peak in the HSQC spectrum based on the average :math:`\tau_c` provided. 
    The function fits the NH region of the 1D spectrum (or of the HSQC reference) to a skew normal distribution - or the integral of the spectrum - to estimate the area of a single peak, which is then used to estimate the intensity of a single peak in the HSQC spectrum. 
    The noise is estimated from a region of the spectrum without peaks. The function prints the estimated SNR per scan for a single peak and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.

    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to load the data and perform the analysis.
        The CO object should have the following attributes:
        
        *   B_0 : float
            The magnetic field strength in Tesla.
        *   S2 : list of float
            The order parameters to use for the calculation of tau_c. In IDP mode, two values should be provide, else only one.
        *   tau : list of float
            The correlation times to use for the calculation of tau_c. In IDP mode, two values should be provide, else only one.
        *   hsqcexpno : str, optional
            The experiment number of the HSQC reference spectrum to use for the SNR estimation.
        *   tractexpno : str, optional
            The experiment number of the TRACT spectrum to use for the SNR estimation. 
        *   T : float, optional
            The temperature in Kelvin, used to recompute tau_c based on temperature.
        *   nres : int, optional
            The number of non-proline residues in the protein. If not provided, it will be estimated from the molecular weight using a rough estimate of 110 Da per residue and multiplying by 0.937 to account for the fact that prolines do not have an amide proton and therefore do not contribute to the signal in the HSQC spectrum.
        

    Returns
    -------
    None
        The function prints the estimated SNR per scan for a single peak, suggests a number of scans for the three experiments, and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.
    """


    if hasattr(CO, 'MW') and CO.MW is not None:
        MW = CO.MW
    else:
        if hasattr(CO, 'lw') and CO.lw is not None:
            print('No value provided for the molecular weight. Estimating it from the linewidth provided. This is a rough estimate, but it should be sufficient for our purposes.')
            MW = (CO.lw*np.pi/7.42 - 0.1674)/0.5998
            CO.add_ref('bermel')
            CO.add_ref('cavanagh')
        elif hasattr(CO, 'nres') and CO.nres is not None:
            MW = CO.nres*0.110/0.937 #average mass of an amino acid is 110 Da, and we multiply by 0.937 to account for the fact that prolines do not have an amide proton and therefore do not contribute to the signal in the HSQC spectrum. This is a rough estimate, but it should be sufficient for our purposes.
        elif hasattr(CO, 'tau') and CO.tau is not None:
            CO.add_ref('cavanagh')
            MW = (hydrodynamics_utils.tau(CO.tau[0]*1e9, 298.15, T_old=CO.T)- 0.1674)/0.5998
        else:   
            raise ValueError('No value provided for the molecular weight or linewidth. Cannot estimate molecular weight.')
    
    if hasattr(CO, 'nres') and CO.nres is not None:
        nres = CO.nres
    else:
        nres = np.rint(MW/0.110)*(0.937) #average mass of an amino acid is 110 Da, and we multiply by 0.937 to account for the fact that prolines do not have an amide proton and therefore do not contribute to the signal in the HSQC spectrum. This is a rough estimate, but it should be sufficient for our purposes.
    
    if CO.options['idp']:
        tau_slow = CO.tau[0]
        tau_int = CO.tau[1]
        S2_slow = CO.S2[0]
        S2_int = CO.S2[1]
        CO.add_ref('rezaei-ghaleh')
        CO.add_ref('salvi')
        CO.add_ref('parigi2000')
        tau_average = S2_slow*tau_slow + S2_int*tau_int
    else:
        tau_average = CO.tau[0]
    if not hasattr(CO, 'lw') or CO.lw is None:
        R2H = tau_average*7.42e9 
        CO.add_ref('bermel')
        amidelinewidth = R2H/(np.pi)    
    else:
        amidelinewidth = CO.lw
    amidelinewidth_p = kz.misc.freq2ppm(amidelinewidth, CO.B_0*kz.sim.gamma['1H'])

    #estimate signal to noise ratio. if HSQC is provided, use it. Else, use the TRACT
    if CO.hsqc is not None:
        path_hsqc = os.path.join(CO.basedir, f'{CO.hsqcexpno}')
        S_hsqc = kz.Spectrum_2D(path_hsqc)
        version = t1t2ne_utils.fs_version(S_hsqc)
        if version == 'topspin3':
            S_hsqc.acqus['GRPDLY'] += 1
        S_hsqc.procs['wf'][-1]['mode'] = 'em'
        S_hsqc.procs['wf'][-1]['lb'] = 1/S_hsqc.acqus['AQ2']
        S_hsqc.procs['zf'][-1] = 2 * S_hsqc.fid.shape[-1]
        S_hsqc.pknl()
        S_hsqc.xf2()
        S_hsqc.projf2(0)
        trace = S_hsqc.Trf2['0.00']
        trace.S = kz.processing.hilbert(trace.r)
        if CO.options['phase']:
            trace.adjph()
        signalregion = kz.fit.get_region(trace.ppm, trace.r, fig_title='Signal region for SNR estimation')
        signal = t1t2ne_utils.extract_regions_from_trace(trace.ppm, trace.r, signalregion)
        noiseregion = kz.fit.get_region(trace.ppm, trace.r, fig_title='Noise region for SNR estimation')        
        noise = t1t2ne_utils.extract_regions_from_trace(trace.ppm, trace.r, noiseregion)
        ntr = S_hsqc.acqus['TD1']
        ns = S_hsqc.ngdic['acqus']['NS']
        x = t1t2ne_utils.extract_regions_from_trace(trace.ppm, trace.ppm, signalregion)
    
    else:
        path = os.path.join(CO.basedir, f'{CO.tract}')


    #load the dataset and check if it's a TRACT experiment
        S = kz.Pseudo_2D(path)
        if not t1t2ne_utils.istract(S):
            raise NameError(f'Experiment {CO.tract} is not a TRACT experiment')
        version = t1t2ne_utils.fs_version(S_hsqc)
        if version == 'topspin3':
            S_hsqc.acqus['GRPDLY'] += 1
        Sa, Sb = split_tract(S)
        if CO.options['phase']:
            Sb.adjph(ref=0)
        signalregion = kz.fit.get_region(Sb.ppm_f2, Sb.rr[0], fig_title='Signal region for SNR estimation')        
        noiseregion = kz.fit.get_region(Sb.ppm_f2, Sb.rr[0], fig_title='Noise region for SNR estimation')
        noise = t1t2ne_utils.extract_regions_from_trace(Sb.ppm_f2, Sb.rr[0], noiseregion)
        signal = t1t2ne_utils.extract_regions_from_trace(Sb.ppm_f2, Sb.rr[0], signalregion)
        ntr = float(input('Enter the number of transients (or press enter to use the default value 128): ').strip() or 128)
        ns = Sb.ngdic['acqus']['NS']
        x = t1t2ne_utils.extract_regions_from_trace(Sb.ppm_f2, Sb.ppm_f2, signalregion)        
    
    CO.add_ref('klassez')    

    result = f_fit.fit_skewnormal(x, signal)
    model, A, a = f_fit.skgaussian_ls(result.params, x, signal, result=True)
    noisestd = kz.anal.noise_std(noise)
    
    print(f'\nEstimated MW of the protein: {MW:.2f} kDa, which corresponds to approximately {nres:.0f} residues')
    #divide the area of the skew normal by the number of residues
    nres = int(input('Enter the number of non-proline residues in the protein (or press enter to estimate it from the molecular weight): ').strip() or nres)
    print(f'Estimated amide proton linewidth: {amidelinewidth:.2f} Hz - {amidelinewidth_p:.2f} ppm.')
    areaofasinglepeak = A/nres
    intensityofasinglepeak = areaofasinglepeak/(amidelinewidth_p*np.sqrt(2*np.pi))
    iosp_perscan = intensityofasinglepeak/(2*np.sqrt(ns))
    noiseinthe2dexperiment = noisestd/np.sqrt(ntr)
    SNR = iosp_perscan/noiseinthe2dexperiment
    print(textcolor(f'\nEstimated SNR per scan for a single peak: {SNR:.2f}\n', 'green'))
    return SNR

def suggest_scans(CO, SNR):
    """
    Suggest a number of scans given the expected reduction in SNR
    
    Parameters    
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to perform the analysis.
        The CO object should have the following attributes:
        
        *   xred : list of float
            The percent residual signal for the longest delay of the experiment, or NOE efficiency as percentage. Accepts a list, computes the optimal number of scans for each value.
    SNR : float
        The estimated signal-to-noise ratio per scan for a single peak.
        
    Returns
    -------
    None
        The function prints the suggested number of scans for each experiment in the list.
    """

    for red in CO.xred:
        if red > 100 or red < 0:
            print(textcolor(f'Invalid value for xred: {red}. It should be between 0 and 100. Skipping this value.', 'red'))
            continue
        ns = int(np.ceil(3/(SNR*red))**2)
        print(textcolor(f'\nFor a percent residual signal of {red*100:.0f}%, the estimated number of scans needed to have an SNR of 3 is {ns}.', 'blue')) 
        print(textcolor(f'The closer multiple of  4 (standard NOE sequence) is {int(np.ceil(ns/4)*4)}', 'blue'))
        print(textcolor(f'The closer multiple of  8 (standard T1, T2 sequence) is {int(np.ceil(ns/8)*8)}', 'blue'))
        print(textcolor(f'The closer multiple of 16 is {int(np.ceil(ns/16)*16)}', 'blue'))       

