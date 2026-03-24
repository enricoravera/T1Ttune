#! /usr/bin/env python3

import numpy as np
import klassez as kz
import argparse
import warnings
from scipy.signal import savgol_filter
from copy import deepcopy
import os
from .textcolor import textcolor
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor

from .base import BaseCommand
from . import f_fit, t1t2ne_utils, fun_hetrelax_models

class TractCmd(BaseCommand):
    SHORT_HELP = "Fit a TRACT experiment to extract the average tau_c of the system"
    DESCRIPTION = '\n'.join([
    'This module contains the functions to analyze a TRACT experiment and extract the average tau_c of the system using the algebraic analysis described in Robson et al. (2021), doi:10.1007/s10858-021-00379-5.',
    'The analysis can be performed using integrals of the spectra or point-by-point fitting, and the spectra can be phased and smoothed if needed. '
    'In the standard operation mode, the tau_c is calculated.'
    'In the IDP mode, the Lipari-Szabo order parameter S2 for the intermediate motion is calculated.'
    'The IDP mode requires input of the molecular weight of the protein to estimate the correlation time of the slow motion, and it assumes a correlation time of 1.6 ns for the intermediate motion and an order parameter of 0.15 for the slow motion, which are typical values for IDPs but might not be accurate for all systems.'
    ,])
    @staticmethod
    def add_arguments(parser: argparse.ArgumentParser) -> None:
        parser.add_argument('--basedir', type=str, help='The base directory where the data is stored.')
        parser.add_argument('--tract', type=str, help='The name of the TRACT experiment to analyze.')
        parser.add_argument('--integrate', '-i', action='store_true', help='Whether to perform the analysis using integrals of the spectra instead of point-by-point fitting.')
        parser.add_argument('--readints', action='store_true', help='Whether to read pre-computed integrals instead of calculating them from the spectra. If --integrate is not used, this option is ignored.')
        parser.add_argument('--selectregion', '-s', action='store_true', help='Whether to select a region of the spectrum for the analysis instead of using a default range.')
        parser.add_argument('--phase', '-p', action='store_true', help='Whether to phase the spectra.')
        parser.add_argument('--S2', nargs='?', help = 'The Lipari-Szabo order parameter S2 to use for the calculation of tau_c. Caution: even if --idp is used, a single value should be provided.')
        parser.add_argument('--idp', action='store_true', help='Whether to use the IDP model to extract the order parameter S2 instead of tau_c.')
        parser.add_argument('--MW', type=float, help='The molecular weight of the protein in kDa, to be used for estimating the tau_slow in the IDP model. If not specified, will raise an error if --idp is used.')
        parser.add_argument('--corr_window_idp', type=float, default=20, help='To estimate the correlation time of the intermediate motion in the IDP model, a correlation time of a peptide of 20 residues is used by default. This parameter allows to change this value if needed.')
        parser.add_argument('--T', type=float, default=298.15, help='The temperature in Kelvin, to be used for estimating the tau_slow in the IDP model. Default is 298 K.')
        parser.add_argument('--slw', type=float, help='The percentage of the spectrum to use for the sliding window smoothing. If not specified, no smoothing is applied.')
        parser.add_argument('--smoothdata', action='store_true', help='Whether to apply a Savitzky-Golay filter to the spectra before the analysis. The window length of the filter is determined by the --slw parameter.')
        parser.add_argument('--smoothrates', action='store_true', help='Whether to apply a Savitzky-Golay filter to the relaxation rates before the analysis. The window length of the filter is determined by the --slw parameter.')
        parser.add_argument('--plot', action='store_true', help='Whether to plot the results of the analysis. If --idp is not used, the tau_c values will be plotted. If --idp is used, the S2 values will be plotted.')   
        parser.add_argument('--nucs', nargs='*', default=['1H', '15N'], help='The nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.')
        parser.add_argument('--r', type=float, default=1.02, help='The length of the 1H-15N bond in Angstroms. Default is 1.02 A.')
        parser.add_argument('--Deltasigma', type=float, default=-160, help='The chemical shift anisotropy of the 15N nucleus in ppm. Default is -160 ppm.')
        parser.add_argument('--theta', type=float, default=17, help='The angle between the 1H-15N bond and the principal axis of the CSA tensor in degrees. Default is 17 degrees.')

    @staticmethod
    def run(args: argparse.Namespace) -> None:
        CO = t1t2ne_utils.Conf_Optns(args, module='tract')
        tract(CO)
        t1t2ne_utils.the_end(CO)
        
def filter_data(y, window_length=5, polyorder=3):
    """Apply a Savitzky-Golay filter to the data y with the specified window length and polynomial order.
    If the data is shorter than the window length, return it unfiltered.
    
    Parameters
    ----------
    y : array_like
        The data to filter.
    window_length : int, optional
        The length of the filter window (i.e., the number of coefficients). Default is 5.
    polyorder : int, optional
        The order of the polynomial used to fit the samples. Default is 3.

    Returns
    -------
    array_like
        The filtered data.
    """
    if y.shape[-1] < window_length:
        return y
    return savgol_filter(y, window_length=window_length, polyorder=polyorder, axis=-1)

def make_plot(CO, xaxis, y_rec, yerr_rec, y_ave, yerr_ave, Ra, Rb, sigma_Ra, sigma_Rb, values, name='tau_c', idx=0):
    """
    Make a plot of the reconstructed tau_c or S2 values as a function of the chemical shift, with error bars and the average value with its error. The plot also includes the expected relaxation rates based on the input parameters, and a reference spectrum of the Sb dataset.
    
    Parameters
    ----------
    CO : Conf_Optns object
        The configuration options object.
    xaxis : array_like
        The x-axis values (chemical shift).
    y_rec : array_like
        The reconstructed tau_c or S2 values.
    yerr_rec : array_like
        The errors of the reconstructed tau_c or S2 values.
    y_ave : float
        The average tau_c or S2 value.
    yerr_ave : float
        The error of the average tau_c or S2 value.
    values : dict
        A dictionary containing the expected relaxation rates and other parameters.
    name : str, optional
        The name of the parameter being plotted ('tau_c' or 'S^2'). Default is 'tau_c'.

    Returns
    -------
    None
    """

    R1 = values['R1']
    R2 = values['R2']
    nOe = values['nOe']
    eta_z = values['eta_z']
    eta_xy = values['eta_xy']
    fig = plt.figure()
    ax = fig.add_subplot(313)
    axr = fig.add_subplot(312, sharex=ax)
    axs = fig.add_subplot(311, sharex=ax)        
    fig.set_size_inches(12, 8)
    axs.plot(CO.Sa.ppm_f2, CO.Sa.rr[idx], label='Sa', color='k')
    axs.plot(CO.Sb.ppm_f2, CO.Sb.rr[idx], label='Sb', color='tab:gray')
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.90)
    axr.errorbar(xaxis, Ra, yerr=sigma_Ra, fmt='.', color='tab:orange', ms=8, ecolor='k', elinewidth=0.4, capsize=2, label=r'$R_a$')
    axr.errorbar(xaxis, Rb, yerr=sigma_Rb, fmt='.', color='tab:green', ms=8, ecolor='k', elinewidth=0.4, capsize=2, label=r'$R_b$')
    #error bars for tau_c
    ax.errorbar(xaxis, y_rec, yerr=yerr_rec, fmt='.', color='tab:blue', ms=8, ecolor='k', elinewidth=0.4, capsize=2, label=r'reconstructed ${}$'.format(name))
    #and average line
    ax.axhline(y=y_ave, color='tab:red', ls='--', zorder=4, label=r'average ${}$ = {:.2e} s $\pm$ {:.2e} s'.format(name, y_ave, yerr_ave))
    #shade the area of the average tau_c +/- its sigma
    ax.fill_between(xaxis, y_ave - yerr_ave, y_ave + yerr_ave, color='tab:red', alpha=0.2)
    #prettify the plot
    kz.misc.pretty_scale(ax, (max(xaxis), min(xaxis)), 'x')
    if name=='tau_c':
        kz.misc.pretty_scale(ax, (1e-9,3e-7), 'y')
    else:
        kz.misc.pretty_scale(ax, (-0.1,1.1), 'y')
    ax.set_xlabel(r'$\delta$ (ppm)')

    if name=='tau_c':
        ax.set_ylabel(r'$\tau_c$ estimate (s)')        
        ax.set_yscale('log')
    else:
        ax.set_ylabel(r'$S^2$ estimate')        
    ax.text(0.05, 0.95, f'expected $T_1$ = {1/R1:.2f} s\nexpected $T_2$ = {1/R2:.2f} s\nexpected hetnOe = {nOe:.2f}\nexpected $\\eta_z$ = {eta_z:.2f}\nexpected $\\eta_{{xy}}$ = {eta_xy:.2f}', transform=ax.transAxes, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    kz.misc.pretty_scale(axr, (0, R2+4*np.abs(eta_xy)), 'y')
    axr.set_ylabel('Relaxation rates (Hz)')
    ax.legend(loc='lower right', bbox_to_anchor=(0.95, 0.905), ncols=2, bbox_transform=fig.transFigure)
    cursor = MultiCursor(fig.canvas, (ax, axs), color='green', linewidth=0.4, horizOn=False)        
    plt.show()

def get_stdev_rate(y, window_length=5):
    """
    Calculate the standard deviation of the data y in a moving window of specified length.
    
    Parameters
    ----------
    y : array_like
        The data to calculate the standard deviation for.
    window_length : int, optional
        The length of the moving window as a percentage of the total data length. Default is 5.

    Returns
    -------
    array_like
        The standard deviation of the data in the moving window.
    """
    
    stdev = np.zeros_like(y)
    count_to = int(len(y) * window_length / 100) // 2
    for j in range(len(y)):
        stdev[j] = np.std(y[max(0, j-count_to):min(y.shape[-1], j+count_to)])
    return stdev

def s2(w_N, c, tau_slow, S2_slow = 0.15, tau_int = 1.6e-9):
    r"""
    This function is applicable for IDPs. The model is based on significant assumptions, and we do not endorse the use of this model for extracting quantitative information.
    The obtained values are only meant to be used to bring the expected relaxation rates in a reasonable ballpark.
    Calculate the Lipari-Szabo order parameter :math:`S^2` for the intermediate motion based on the work of Razaei-Ghaleh et al. (2018), `doi:10.1002/anie.201808172`_. 

    .. _doi:10.1002/anie.201808172: https://doi.org/10.1002/anie.201808172
    
    The slow motion is assumed to correspond to the tumbling of a protein of the same size of the one under study, and the order parameter of the slow motion is assumed to be 0.15.
    The intermediate motion is assumed to have a correlation time of a peptide of 20 residues (1.6 ns).
    
    Parameters
    ----------
    w_N : float
        The Larmor frequency of the nucleus of interest (in radians per second).
    c : float
        half of the difference of Ra and Rb, to which :math:`\eta_{xy}/(d*c*B_0*P_2(\cos(\theta)))` should be equal
    tau_slow : float
        The correlation time of the slow motion (in seconds).
    S2_slow : float, optional
        The order parameter of the slow motion. Default is 0.15.
    tau_int : float, optional
        The correlation time of the intermediate motion (in seconds). Default is 1.6e-9.
        
    Returns
    -------
    float
        The Lipari-Szabo order parameter :math:`S^2` for the intermediate motion.

    """
    
    J_0_slow = fun_hetrelax_models.J_omega_tau_iso(0, tau_slow)
    J_nuc2_slow = fun_hetrelax_models.J_omega_tau_iso(w_N, tau_slow)
    J_0_intermediate = fun_hetrelax_models.J_omega_tau_iso(0, tau_int)
    J_nuc2_intermediate = fun_hetrelax_models.J_omega_tau_iso(w_N, tau_int)
    J_0_fast = fun_hetrelax_models.J_omega_tau_iso(0, 1e-11)
    J_nuc2_fast = fun_hetrelax_models.J_omega_tau_iso(w_N, 1e-11)
    
    num = (4) * (S2_slow * J_0_slow + (1 - S2_slow) * J_0_fast) + (3) * (S2_slow * J_nuc2_slow + (1 - S2_slow) * J_nuc2_fast)
    denom = (4) * (1-S2_slow) * (J_0_intermediate - J_0_fast) + (3) * (1-S2_slow) * (J_nuc2_intermediate - J_nuc2_fast)

    S2 = (np.abs(c) - num) / (denom)    

    return S2

def tc(w_N, c, O):
    r"""
    Calculate the optimal :math:`\tau_c` based on the equation 15 from Robson et al. (2021), `doi:10.1007/s10858-021-00379-5`_.
    
    .. _doi:10.1007/s10858-021-00379-5: https://doi.org/10.1007/s10858-021-00379-5
    
    Parameters
    ----------
    w_N : float
        The Larmor frequency of the nucleus of interest (in radians per second).
    c : float
        half of the difference of Ra and Rb, to which :math:`\eta_{xy}/(d*c*B_0*P_2(\cos(\theta)))` should be equal    
    O : float
        Lipari-Szabo order parameter :math:`S^2`.
        
    Returns
    -------
    float
        The closed form solution for :math:`\tau_c` (in seconds).
    """
    
    t1 = (5*c)/(24*O)
    D = 24*O*w_N**2
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        Rq = np.sqrt(21952*(O**6)*(w_N**6) - 3025*(O**4)*(c**2)*(w_N**8) + 625*(O**2)*(c**4)*(w_N**10))
        Rc = (1800*(O**2)*c*(w_N**4) + 125*(c**3)*(w_N**6) + 24*np.sqrt(3)*Rq)**(1/3)
        t2 = (336*(O**2)*(w_N**2) - 25*(c**2)*(w_N**4)) / (D * Rc)
        t3 = Rc / D
    
    #t1 = (5*c)/24 
    #t2 = (336*(w_N**2) - 25*(c**2)*(w_N**4)) / (24*(w_N**2) * (1800*c*(w_N**4) + 125*(c**3)*(w_N**6) + 24*np.sqrt(3)*np.sqrt(21952*(w_N**6) - 3025*(c**2)*(w_N**8) + 625*(c**4)*(w_N**10)))**(1/3))
    #t3 = (1800*c*(w_N**4) + 125*(c**3)*(w_N**6) + 24*np.sqrt(3)*np.sqrt(21952*(w_N**6) - 3025*(c**2)*(w_N**8) + 625*(c**4)*(w_N**10)))**(1/3)/(24*w_N**2) 
    
    return t1 - t2 + t3    

def tract_fit_Ra_Rb(CO):
    r"""
    Fit the TRACT experiment to extract the relaxation rates :math:`R_a` and :math:`R_b` of the two components of the TRACT experiment using the algebraic analysis described in Robson et al. (2021), `doi:10.1007/s10858-021-00379-5`_.
    
    .. _doi:10.1007/s10858-021-00379-5: https://doi.org/10.1007/s10858-021-00379-5
    
    The monoexponential fit is performed using :func:`fit_exponential` from `TRAGICO code`_.
    
    .. _TRAGICO code: https://github.com/letiziafiorucci/tragico  
    
    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to load the data and perform the analysis.
        CO must contain the following options:
        
        *   integrate : bool
            Whether to perform the analysis using integrals of the spectra instead of point-by-point fitting. Default is False.
        *   selectregion : bool
            Whether to select a region of the spectrum for the analysis instead of using a default range. Default is False.
        *   phase : bool
            Whether to phase the spectra. Default is False.
        *   basedir : str
            The base directory where the data is stored.
        *   tract : str
            The name of the TRACT experiment to analyze.

        CO must also contain the following attributes:

        *   r : float
            The length of the 1H-15N bond in meters. Default is 1.02e-10 m.
        *   Deltasigma : float
            The chemical shift anisotropy of the 15N nucleus. Default is -160 ppm.
        *   theta : float
            The angle between the 1H-15N bond and the principal axis of the CSA tensor in radians. Default is 17 degrees converted to radians.
        *   nucs : list of str
            The list of nuclei to use for the calculation of the relaxation rates. Default is ['1H', '15N'].  
    Returns
    -------
    tuple
        A tuple containing the x axis for plotting, the arrays of :math:`R_a` and :math:`R_b` values for each point in the spectrum, and the corresponding uncertainties.  
        
    """
    
    CO.add_ref('lee')
    CO.add_ref('robson')
    CO.add_ref('Fiorucci')
    CO.add_ref('klassez')

    path = os.path.join(CO.basedir, f'{CO.tract}')
    # find vdlist in path

    #   Connect the name of the file to the path
    S = kz.Pseudo_2D(path)
    if not t1t2ne_utils.istract(S):
        raise NameError(f'Experiment {CO.tract} is not a TRACT experiment')
    try:
        file = os.listdir(os.path.join(path, 'lists', 'vd'))[0]
        vdlistpath = os.path.join(path, 'lists', 'vd', file)
    except:
        vdlistpath = os.path.join(path, 'vdlist')
    #   Print a notification
    print(f'Found {vdlistpath} to be imported as VDLIST')
    #   Actual loading and storage in an attribute
    vdlist = t1t2ne_utils.in_vdlist(vdlistpath)
    print(f'vdlist loaded: {vdlist}')
    
    idx = np.argmin(np.abs(vdlist))
    
    integrate = CO.options['integrate']
    readints = CO.options['readints']
    selectregion = CO.options['selectregion']
    phase = CO.options['phase']
    smooth_data = CO.options['smoothdata']
    smooth_rates = CO.options['smoothrates']
    if hasattr(CO.options, 'slw') and CO.options['slw'] is not None:
        slw = CO.options['slw']

    #load the dataset and check if it's a TRACT experiment

    CO.B0 = S.acqus['SFO1'] / kz.sim.gamma[CO.nucs[0]]
    print(f'Magnetic field strength: {CO.B0:.2f} T')
    #splitcomb-like operation to separate the two interleaved datasets (TROSY and ANTITROSY)
    version = t1t2ne_utils.fs_version(S)
    if version == 'topspin4':
        S.fid = np.reshape(S.fid.flatten(), (2*S.fid.shape[0], -1))
    else:
        pass
    Sa = deepcopy(S)
    Sb = deepcopy(S)
    Sb.fid = S.fid[::2]
    Sa.fid = S.fid[1::2]
    Sa.acqus['TD1'] = Sa.fid.shape[0]
    Sb.acqus['TD1'] = Sb.fid.shape[0]
    if smooth_data or smooth_rates:
        slw_windowlength = int(S.acqus['TD'] * slw[0] / 100)

    #phase the two datasets 
    for s in [Sa, Sb]:
        # Window function: exponential with 0.5 Hz linebroadening
        s.procs['wf']['mode'] = 'em'
        s.procs['wf']['lb'] = 1/s.acqus['AQ']
        if s==Sa:
            print(f'Line broadening applied: {s.procs["wf"]["lb"]:.2g} Hz')
        # Zerofill to twice the size
        s.procs['zf'] = 2 * s.fid.shape[-1]
        # Apply and do FT
        if phase:
            s.procs['p0'] = 0
            s.procs['p1'] = 0
        s.process()

        # Remove digital filter
        s.pknl()

        # Phase the spectrum
        if phase:

            if s == Sa:
                s.adjph()
            else:
                s.adjph(p0=Sa.procs['p0'], p1=Sa.procs['p1'], pv=Sa.procs['pv'], update=False)
        # Qfil the solvent peak
        s.acqus['FnMODE'] = 'No'
        s.abs_v2()
        if s==Sa:
            s.qfil()
        else:
            s.qfil(u=Sa.procs['qfil']['u'], s=Sa.procs['qfil']['s'])
        if smooth_data:
            s.rr = filter_data(s.rr, window_length=slw_windowlength)


    Rb = []
    Ra = []
    sigma_Rb = []
    sigma_Ra = []

    if integrate:

    # Analysis using integrals of the two spectra
        filenames = 's_b', 's_a'

    #       Integrals
        for k, (filename, s) in enumerate(zip(filenames, [Sa, Sb])):
            if k == 0:    # => TROSY component
                if not readints:       # compute the integrals  
                    s.integrate(filename=filename)
                else:       # read an integrals file
                    s.read_integrals(filename=f'{filename}.igrl')
            else:       # => ANTITROSY component
                # Get the integration regions from the TROSY spectrum
                limits = kz.misc.key_to_limits(list(Sa.integrals.keys()))
                if not readints:       # integrate in the same regions of the TROSY
                    s.integrate(filename=filename, lims=limits)
                else:       # read an integrals file
                    s.read_integrals(filename=f'{filename}.igrl')
        for key in Sb.integrals.keys():
            #fit the two datasets to get the relaxation rates Rb and Ra, which are the rates of the two components of the TRACT experiment.
            
            Rb.append(f_fit.fit_exponential(vdlist, Sb.integrals[key], multi=1).params['k'].value)
            Ra.append(f_fit.fit_exponential(vdlist, Sa.integrals[key], multi=1).params['k'].value)
            sigma_Rb.append(f_fit.fit_exponential(vdlist, Sb.integrals[key], multi=1).params['k'].stderr)
            sigma_Ra.append(f_fit.fit_exponential(vdlist, Sa.integrals[key], multi=1).params['k'].stderr) 

        #print(limits)
        if len(np.asarray(limits).shape) == 1:
            xaxis = [np.mean(np.asarray(limits))]
        else:   
            xaxis = np.mean(np.asarray(limits), axis=1)

    else:
        xaxis = []
        if selectregion:
            signalregion = kz.fit.get_region(Sa.ppm_f2, Sa.rr[0])
        else:
            if CO.options['idp']:
                signalregion = [[8.6, 7.4]]
            else:
                signalregion = [[10, 7]]
        # Sb.strip(signalregion)
        # Sa.strip(signalregion)
        for sreg in signalregion:
            start, end = [kz.misc.ppmfind(Sa.ppm_f2, w)[0] for w in sreg]
            for i in range(start, end):
                Rb.append(f_fit.fit_exponential(vdlist, Sb.rr[:, i], multi=1).params['k'].value)
                Ra.append(f_fit.fit_exponential(vdlist, Sa.rr[:, i], multi=1).params['k'].value)
                sigma_Rb.append(f_fit.fit_exponential(vdlist, Sb.rr[:, i], multi=1).params['k'].stderr)
                sigma_Ra.append(f_fit.fit_exponential(vdlist, Sa.rr[:, i], multi=1).params['k'].stderr)
                xaxis.append(Sb.ppm_f2[i])
    Rb = [np.exp(r) for r in Rb]
    Ra = [np.exp(r) for r in Ra]
    sigma_Rb = [r * s for r, s in zip(Rb, sigma_Rb)]
    sigma_Ra = [r * s for r, s in zip(Ra, sigma_Ra)]

    Rb = np.array(Rb)
    Ra = np.array(Ra)
    if smooth_rates:
        Rb = filter_data(Rb, window_length=slw_windowlength)
        Ra = filter_data(Ra, window_length=slw_windowlength)
    sigma_Rb = np.array(sigma_Rb)
    sigma_Ra = np.array(sigma_Ra)
    if smooth_rates:
        sigma_Rb = get_stdev_rate(Rb, window_length=slw_windowlength)
        sigma_Ra = get_stdev_rate(Ra, window_length=slw_windowlength)
    CO.Sb = Sb
    CO.Sa = Sa
    CO.xaxis = xaxis
    #CO.signalregion = signalregion
    
    return xaxis, Ra, Rb, sigma_Ra, sigma_Rb, S.acqus['B0'], idx

def tract_compute_tau(B0, Ra, Rb, sigma_Ra, sigma_Rb, S2=0.9, r=1.02e-10, Deltasigma=-160, theta=17*np.pi/180, nuc1='1H', nuc2='15N'):
    """
    Compute the tau_c values from the relaxation rates Ra and Rb obtained from the TRACT experiment using the algebraic analysis described in Robson et al. (2021), `doi:10.1007/s10858-021-00379-5`_.
    
    .. _doi:10.1007/s10858-021-00379-5: https://doi.org/10.1007/s10858-021-00379-5
    
    Parameters
    ----------

    w_N : float
        The Larmor frequency of the nucleus of interest (in radians per second).
    Ra : array_like
        The relaxation rates of the ANTITROSY component of the TRACT experiment.
    Rb : array_like
        The relaxation rates of the TROSY component of the TRACT experiment.
    sigma_Ra : array_like
        The standard deviations of the relaxation rates of the ANTITROSY component.
    sigma_Rb : array_like
        The standard deviations of the relaxation rates of the TROSY component.
    r : float, optional
        The length of the 1H-15N bond in meters. Default is 1.02e-10 m.
    Deltasigma : float, optional
        The chemical shift anisotropy of the 15N nucleus. Default is -160 ppm.
    theta : float, optional
        The angle between the 1H-15N bond and the principal axis of the CSA tensor in radians. Default is 17 degrees converted to radians.
    nuc1 : str, optional
        The name of the first nucleus involved in the dipolar interaction (default is '1H').
    nuc2 : str, optional
        The name of the second nucleus involved in the dipolar interaction (default is '15N').
        
    Returns
    -------
    tuple
        A tuple containing the array of tau_c values for each point in the spectrum and the corresponding uncertainties, the average tau_c and its uncertainty.
    """
    # equations (5) and (6) of https://pmc.ncbi.nlm.nih.gov/articles/PMC8627365/
    #sqrt2 for both p and dN is not included and is not needed since it cancels out in the final equation
    p = fun_hetrelax_models.d(r=r, nuc1=nuc1, nuc2=nuc2) # DD 1H-15N bond in rad/s
    # equation (6)
    dN = fun_hetrelax_models.c(Deltasigma=Deltasigma, nuc=nuc2) * B0        # 15N CSA
    # equation (7)
    w_N = B0 * kz.sim.gamma[nuc2] * 2 * np.pi * 1e6                # 15N frequency (radians/s)
        # from equivalence between equation (8) RHS and equation (9) RHS

    c = (Rb - Ra)/(dN*p*(3*np.cos(theta)**2-1))

    #compute the uncertainty on c by propagating the uncertainties on Rb and Ra, assuming they are independent

    sigma_eta = 0.5 * np.sqrt(sigma_Rb**2 + sigma_Ra**2)
    dc_deta = 1 / (dN*p*(3*np.cos(theta)**2-1))
    sigma_c = np.abs(dc_deta) * sigma_eta

    #compute tau_c and its uncertainty by propagating the uncertainty on c using numerical differentiation, since the equation is complicated and it's easier to do it numerically than analytically.
    tau_c = tc(w_N, c, S2)
 
    dec = 1e-3*sigma_c
    detau_dec = (tc(w_N, c+dec, S2) - tc(w_N, c-dec, S2))/(2*dec)
    sigma_tau_c = np.abs(detau_dec) * sigma_c
    print('Dipolar coupling constant p = {:.2e} rad/s'.format(p))
    print('CSA constant dN = {:.2e} rad/s'.format(dN))
    print('Larmor frequency w_N = {:.2e} rad/s'.format(w_N))
    print('angle theta = {:.2f} degrees'.format(theta*180/np.pi))
    Rb_avg = np.mean(Rb)
    Ra_avg = np.mean(Ra)
    sigma_Rb_avg = np.std(Rb)
    sigma_Ra_avg = np.std(Ra)
    c_avg = np.average(c, weights=c**2/sigma_c**2)  #(Rb_avg - Ra_avg)/(dN*p*(3*np.cos(theta)**2-1))
    sigma_eta_avg = 0.5 * np.sqrt(sigma_Rb_avg**2 + sigma_Ra_avg**2)
    sigma_c_avg = np.abs(dc_deta) * sigma_eta_avg
    dec_avg = 1e-3*sigma_c_avg
    tau_average = tc(w_N, c_avg, S2)

    detau_dec_avg = (tc(w_N, c_avg+dec_avg, S2) - tc(w_N, c_avg-dec_avg, S2))/(2*dec_avg)
    sigma_tau_average = np.abs(detau_dec_avg) * sigma_c_avg
    
    print(textcolor(f'average Rb = {Rb_avg:.2f} s^-1 ± {sigma_Rb_avg:.2f} s^-1', 'blue', bold=False))
    print(textcolor(f'average Ra = {Ra_avg:.2f} s^-1 ± {sigma_Ra_avg:.2f} s^-1', 'blue', bold=False))
    #compute the average tau_c and its uncertainty using weighted average
    print(textcolor(f'average c = {c_avg:.2e} ± {sigma_c_avg:.2e}', 'blue', bold=False))
    print(textcolor(f'average tau_c = {tau_average:.2e} s ± {sigma_tau_average:.2e} s', 'blue', bold=False))

    return tau_c, sigma_tau_c, tau_average, sigma_tau_average

def tract_compute_s2(B0, Ra, Rb, sigma_Ra, sigma_Rb, tau_slow, S2_slow = 0.15, tau_int = 1.6e-9, r=1.02e-10, Deltasigma=-160, theta=17*np.pi/180, nuc1='1H', nuc2='15N'):
    """
    Compute the Lipari-Szabo order parameter :math:`S^2` for the intermediate motion based on the work of Razaei-Ghaleh et al. (2018), `doi:10.1002/anie.201808172`_.
    
    .. _doi:10.1002/anie.201808172: https://doi.org/10.1002/anie.201808172
    
    The slow motion is assumed to correspond to the tumbling of a protein of the same size of the one under study, and the order parameter of the slow motion is assumed to be 0.15.
    The intermediate motion is assumed to have a correlation time of a peptide of 20 residues (1.6 ns).
    
    Parameters
    ----------
    
    w_N : float
        The Larmor frequency of the nucleus of interest (in radians per second).
    Ra : array_like
        The relaxation rates of the ANTITROSY component of the TRACT experiment.
    Rb : array_like
        The relaxation rates of the TROSY component of the TRACT experiment.
    sigma_Ra : array_like
        The standard deviations of the relaxation rates of the ANTITROSY component.
    sigma_Rb : array_like
        The standard deviations of the relaxation rates of the TROSY component.
    tau_slow : float
        The correlation time of the slow motion (in seconds) estimated from MW.
    S2 : list of a single float, optional
        The order parameter of the slow motion. Default is 0.15.
    tau_int : float, optional
        The correlation time of the intermediate motion (in seconds). Default is 1.6e-9.
    r : float, optional
        The length of the 1H-15N bond in meters. Default is 1.02e-10 m.
    Deltasigma : float, optional
        The chemical shift anisotropy of the 15N nucleus. Default is -160 ppm.
    theta : float, optional
        The angle between the 1H-15N bond and the principal axis of the CSA tensor in radians. Default is 17 degrees converted to radians.
    nuc1 : str, optional
        The name of the first nucleus involved in the dipolar interaction (default is '1H').
    nuc2 : str, optional
        The name of the second nucleus involved in the dipolar interaction (default is '15N').
    
    Returns
    -------
    tuple
        The array of the Lipari-Szabo order parameter :math:`S^2` for the intermediate motion and its uncertainty, the average :math:`S^2` and its uncertainty.
    """
    # equations (5) and (6) of https://pmc.ncbi.nlm.nih.gov/articles/PMC8627365/
    #sqrt2 for both p and dN is not included and is not needed since it cancels out in the final equation
    p = fun_hetrelax_models.d(r=r, nuc1=nuc1, nuc2=nuc2) # DD 1H-15N bond in rad/s
    # equation (6)
    dN = fun_hetrelax_models.c(Deltasigma=Deltasigma, nuc=nuc2) * B0        # 15N CSA
    # equation (7)
    w_N = B0 * kz.sim.gamma[nuc2] * 2 * np.pi * 1e6                # 15N frequency (radians/s)
        # from equivalence between equation (8) RHS and equation (9) RHS
    c = (Rb - Ra)/(dN*p*(3*np.cos(theta)**2-1))
    #compute the uncertainty on c by propagating the uncertainties on Rb and Ra, assuming they are independent
    sigma_eta = 0.5 * np.sqrt(sigma_Rb**2 + sigma_Ra**2)
    dc_deta = 1 / (dN*p*(3*np.cos(theta)**2-1))
    sigma_c = np.abs(dc_deta) * sigma_eta

    S2 = s2(w_N, c, tau_slow, S2_slow, tau_int)
    #compute the uncertainty on S2 by propagating the uncertainty on c using numerical differentiation, since the equation is complicated and it's easier to do it numerically than analytically.
    dec = 1e-3*sigma_c
    dS2_dec = (s2(w_N, c+dec, tau_slow, S2_slow, tau_int) - s2(w_N, c-dec, tau_slow, S2_slow, tau_int))/(2*dec)
    sigma_S2 = np.abs(dS2_dec) * sigma_c
    print('Dipolar coupling constant p = {:.2e} rad/s'.format(p))
    print('CSA constant dN = {:.2e} rad/s'.format(dN))
    print('Larmor frequency w_N = {:.2e} rad/s'.format(w_N))
    print('angle theta = {:.2f} degrees'.format(theta*180/np.pi))
    
    Rb_avg = np.mean(Rb)
    Ra_avg = np.mean(Ra)
    sigma_Rb_avg = np.std(Rb)
    sigma_Ra_avg = np.std(Ra)
    c_avg = np.average(c, weights=c**2/sigma_c**2)
    sigma_eta_avg = 0.5 * np.sqrt(sigma_Rb_avg**2 + sigma_Ra_avg**2)
    sigma_c_avg = np.abs(dc_deta) * sigma_eta_avg
    dec_avg = 1e-3*sigma_c_avg
    S2_average = s2(w_N, c_avg, tau_slow, S2_slow, tau_int)
    dS2_dec_avg = (s2(w_N, c_avg+dec_avg, tau_slow, S2_slow, tau_int) - s2(w_N, c_avg-dec_avg, tau_slow, S2_slow, tau_int))/(2*dec_avg)
    sigma_S2_average = np.abs(dS2_dec_avg) * sigma_c_avg
    
    print(textcolor(f'average Rb = {np.mean(Rb):.2f} s^-1 ± {np.std(Rb):.2f} s^-1', 'blue', bold=False))
    print(textcolor(f'average Ra = {np.mean(Ra):.2f} s^-1 ± {np.std(Ra):.2f} s^-1', 'blue', bold=False))

    print(textcolor(f'average S2 = {np.mean(S2):.2f} ± {np.mean(sigma_S2):.2f}', 'blue', bold=False))
    return S2, sigma_S2, np.mean(S2), np.mean(sigma_S2)

def tract(CO):

    xaxis, Ra, Rb, sigma_Ra, sigma_Rb, B0_exp, idx = tract_fit_Ra_Rb(CO)
    if CO.options['idp']:
        S2, sigma_S2, S2_average, sigma_S2_average = tract_compute_s2(B0_exp, Ra, Rb, sigma_Ra, sigma_Rb, tau_slow=CO.tau[0], S2_slow=CO.S2[0], tau_int=CO.tau[1], r=CO.r, Deltasigma=CO.Deltasigma, theta=CO.theta, nuc1=CO.nucs[0], nuc2=CO.nucs[1])
        CO.S2[1]=S2_average
    else:
        tau_c, sigma_tau_c, tau_average, sigma_tau_average = tract_compute_tau(B0_exp, Ra, Rb, sigma_Ra, sigma_Rb, S2=CO.S2[0], r=CO.r, Deltasigma=CO.Deltasigma, theta=CO.theta, nuc1=CO.nucs[0], nuc2=CO.nucs[1])    
        CO.tau=[tau_average, 1e-11]
    CO.add_ref('fushman')
    R1, R2, nOe = fun_hetrelax_models.R1R2nOe(B0_exp, r=CO.r, Deltasigma=CO.Deltasigma, nuc1=CO.nucs[0], nuc2=CO.nucs[1], func=fun_hetrelax_models.LS_iso, f_args=(CO.S2, CO.tau))
    CO.add_ref('salvi')
    eta_z, eta_xy = fun_hetrelax_models.eta_z_eta_xy(B0_exp, r=CO.r, Deltasigma=CO.Deltasigma, theta=CO.theta, nuc1=CO.nucs[0], nuc2=CO.nucs[1], f_args=(CO.S2, CO.tau))
    print(f'expected T1 = {1/R1:.3f} s\nexpected T2 = {1/R2:.3f} s\nexpected hetnOe = {nOe:.3f} \nexpected eta_z = {eta_z:.3f}\nexpected eta_xy = {eta_xy:.3f}')
    print('\n')

    if CO.options['plot']:
        values = {'R1': R1, 'R2': R2, 'nOe': nOe, 'eta_z': eta_z, 'eta_xy': eta_xy} 
        if not CO.options['idp']:
            y = tau_c
            yerr = sigma_tau_c
            y_ave = tau_average
            yerr_ave = sigma_tau_average
            name = 'tau_c'
        else:
            y = S2
            yerr = sigma_S2
            y_ave = S2_average
            yerr_ave = sigma_S2_average
            name = 'S^2'
        make_plot(CO, xaxis, y, yerr, y_ave, yerr_ave, Ra, Rb, sigma_Ra, sigma_Rb, values, name=name, idx = idx)
    
    print(textcolor('Use this command to create the lists for your system', 'green'))
    print(f't1t2ne makelists --tau {CO.tau[0]*1e9:.2e} {CO.tau[1]*1e9:.2e} --S2 {CO.S2[0]:.2f} {CO.S2[1]:.2f} --idp' if CO.options['idp'] else f't1t2ne makelists --tau {CO.tau[0]*1e9:.2e} --S2 {CO.S2[0]:.2f}')     
