#! /usr/bin/env python3

import numpy as np

from .base import BaseCommand
from .textcolor import textcolor
from . import t1ttune_utils, fun_hetrelax_models
import random

class SetupTractCmd(BaseCommand):
    SHORT_HELP = "Setup the vdlist for a TRACT experiment"
    DESCRIPTION = '\n'.join([
    'With this module, the vdlist for a TRACT experiment is set up.' 
    'In the standard operation mode, the tau_c or MW are used to compute the optimal vdlist.'
    'S2 is assumed to be 0.9 for the standard operation mode.'
    'The IDP mode requires input of the molecular weight of the protein to estimate the correlation time of the slow motion.',
    'Optionally, the length of the correlation window can be provided as an estimate for the intermediate correlation time.',
    'If not provided, it assumes a correlation time of 1.6 ns for the intermediate motion.'
    'Two S2 values are expected in the IDP mode. If not provided, they are assumed to be 0.15 for the slow motion and 0.35 for the intermediate motion.'
    ,])
    
    @staticmethod
    def add_arguments(parser):
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
        parser.add_argument('--Larmor', nargs=1, type=float, help='The Larmor frequency in MHz, alternative to --B0.')
        parser.add_argument('--nucs', nargs='*', default=['1H', '15N'], help='The nuclei to use for the calculation of the relaxation rates. Default is 1H and 15N.')
        parser.add_argument('--large', action='store_true', help='Whether to create the lists for the "large" sequence, which is optimized for short T2 times. If True, the d21 value is set to 450 us and only 8 cycles per CPMG block are used instead of 16. Default is False.')
        parser.add_argument('--small', action='store_true', help='Whether to create the lists for the ".idp" sequence, which is optimized for long T2 times. If True, the d21 value is set to 600 us. Default is False.')
        parser.add_argument('--randomize', action='store_true', help='Whether to randomize the order of the values in the lists. Default is False.')

    @staticmethod
    def run(args):
        CO = t1ttune_utils.Conf_Optns(args, module='setuptract')
 
        suggest_tract_vdlist(CO)
        t1ttune_utils.the_end(CO)    
        exit()        

def suggest_tract_vdlist(CO):
    r"""
    Suggests a vdlist for the TRACT experiment based on the computed average :math:`\tau_c` from the MW. The longest delay of the vdlist is set to 2 times the Rb relaxation time, which corresponds to a reduction of the signal intensity of about 86%. The vdlist is logarithmically spaced from 20 us to the longest delay. The user can choose the number of increments in the vdlist. The function prints the suggested vdlist in Bruker-readable format and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.

    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the necessary information to load the data and perform the analysis.
        It is expected that the ``CO`` object has the following attributes:
        - ``B_0``: the magnetic field strength in Tesla
        - ``tau``: the average correlation time in seconds, calculated from the MW or provided by the user - list of 3 values for the IDP model, 2 values for the non-IDP model
        - ``S2``: the order parameter S2 to use for the calculation - list of 2 values for the IDP model, 1 value for the non-IDP model
        - ``nucs``: the nuclei to use for the calculation of the relaxation rates. Default is ['1H', '15N'].  
        - ``r``: the length of the 1H-15N bond in meters. Default is 1.02e-10 m.
        - ``Deltasigma``: the chemical shift anisotropy of the 15N nucleus. Default is -160 ppm.
        - ``theta``: the angle between the 1H-15N bond and the principal axis of the CSA tensor in radians. Default is 17 degrees converted to radians.
        - ``nT``: the number of increments in the suggested vdlist. Default is 8.
    
    Returns
    -------
    None
         The function prints the suggested vdlist in Bruker-readable format and updates the ``references`` attribute of the ``CO`` object with the references used for the calculations.

    """
    if hasattr(CO, 'nT'):
        nT = CO.nT[0] if isinstance(CO.nT, list) else CO.nT
    else:
        nT = 8
    CO.add_ref('lee')
    if CO.options['idp']:
        CO.add_ref('rezaei-ghaleh')
    CO.add_ref('fushman')
    R1, R2, nOe = fun_hetrelax_models.R1R2nOe(CO.B_0, r=CO.r, nuc1=CO.nucs[0], nuc2=CO.nucs[1], Deltasigma=CO.Deltasigma, func=fun_hetrelax_models.LS_iso, f_args=(CO.S2, CO.tau))
    CO.add_ref('salvi')
    eta_z, eta_xy = fun_hetrelax_models.eta_z_eta_xy(CO.B_0, r=CO.r, nuc1=CO.nucs[0], nuc2=CO.nucs[1], Deltasigma=CO.Deltasigma, theta=CO.theta, func=fun_hetrelax_models.LS_iso, f_args=(CO.S2, CO.tau))
    Rb = R2 + eta_xy
    Ra = R2 - eta_xy
    vdlist_TRACT = np.geomspace(2e-5, 2/Rb, num=nT) #geometrically spaced list from 20us to 2*tau_average
    if CO.options['randomize']:
        random.shuffle(vdlist_TRACT)
    print(textcolor('\nSuggested vdlist for TRACT experiment:', 'blue'))
    print('-'*38)
    t1ttune_utils.out_vdlist(vdlist_TRACT)
    
