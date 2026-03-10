#! /usr/bin/env python3

import os
import socket
import configparser
import numpy as np

from . import hydrodynamics_utils, f_findfs


def config_exists():
    """
    Checks if the config file exists in the curdir folder of topspin for any of the users in the nmrsuperuser list. 
    If it exists, it returns the path to the config file, otherwise it returns None.
    
    Returns
    -------
    str or None
        The path to the config file if it exists, otherwise None.
    """
    findfs_info = f_findfs.find_topspin()
    if not findfs_info["found"]:
        print("TopSpin does not appear to be installed on this workstation.")
        return None
    if not findfs_info["spectrometer"]:
        print("TopSpin does not appear to be installed on a spectrometer workstation.")
        return None
    fspath = findfs_info["install_path"]
    nmrsupath = os.path.join(fspath, 'conf', 'nmrsuperuser')
    with open(nmrsupath) as file:
        nmrsulist = [line.rstrip() for line in file]
    for user in nmrsulist:
        pathtocurdir = os.path.join(fspath, 'prog', 'curdir', user, 'tract_analysis_config.ini')
        if os.path.exists(pathtocurdir):
            return pathtocurdir
    return None

def load_config():
    """
    Loads the config file from the curdir folder of topspin for the current user and returns a configparser object. 
    If the config file does not exist, it raises a FileNotFoundError.
    
    Returns
    -------
    configparser.ConfigParser
        The loaded config object.
    or
    None
    """
    pathtocurdir = config_exists()
    if pathtocurdir is None:
        print('Config file not found in the curdir folder of topspin for any of the users in the nmrsuperuser list.')
        return None
    else:
        config_p = configparser.ConfigParser()
        config_p.read(pathtocurdir)
        return config_p


def extract_regions_from_trace(scale, trace, region):
    """
    Extract regions of the spectrum and return a single array
    
    Parameters
    ----------
    scale : array_like
        The scale corresponding to the trace (e.g., ppm values).
    trace : array_like
        The spectrum to extract the regions from.
    region : list of tuples
        The regions to extract, each tuple contains the start and end ppm values.

    Returns
    -------
    array_like
        The extracted regions concatenated into a single array.
    """
    for i, reg in enumerate(region):
        if i==0:
            extracted = trace[(scale < max(reg)) & (scale > min(reg))]
        else:
            extracted = np.concatenate((extracted, trace[(scale < max(reg)) & (scale > min(reg))]))
    return extracted

def f4(value):
    """
    Format floats for flopspin. If the value is smaller than 5e-4, it is formatted in micro (u), if it is smaller than 5e-2, it is formatted in milli (m), otherwise it is formatted in the original unit. The value is always formatted with 5 decimal places.

    Parameters
    ----------
    value : float
        The value to be formatted.

    Returns
    -------
    str
        The formatted string.
    """
    if value < 5e-4:
        return f'{value*1e6:.2f}u'
    if value < 5e-2:
        return f'{value*1e3:.2f}m'
    return f'{value:.5f}'

def istract(S):
    """
    Check if the experiment is a TRACT experiment by looking at the name of the sequence 
    (PULPROG parameter in the acqus dictionary of the spectrum).
    
    Parameters    
    ----------
    S : kz.Pseudo_2D
        The spectrum to check.
        
    Returns
    -------
    bool
        True if the experiment is a TRACT experiment, False otherwise.
    """
    
    if 'tract' not in S.ngdic['acqus']['PULPROG'].lower():
        return False
    else:
        return True

def out_vdlist(vdlist):
    """
    Print vdlist in Bruker-readable format.

    Parameters
    ----------
    vdlist : list of float
        The list of values to be printed.
    """
    
    for value in vdlist:
        print(f4(value))

def in_vdlist(vdlist_file):
    """
    Load vdlist from a Bruker-readable file.

    Parameters
    ----------
    vdlist_file : str
        The path to the file containing the vdlist.

    Returns
    -------
    numpy.ndarray
        The array of values read from the file.
    """
    vdlist = []
    with open(vdlist_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.endswith('u'):
                vdlist.append(float(line[:-1]) * 1e-6)
            elif line.endswith('m'):
                vdlist.append(float(line[:-1]) * 1e-3)
            else:
                vdlist.append(float(line))
    return np.array(vdlist)

def splashscreen(module = None):
    """
    Display the splash screen with the name of the software and the authors. If a module name is provided, it is also displayed in the splash screen.
    
    Parameters
    ----------
    module : str, optional
        The name of the module to be displayed in the splash screen. Default is None.
    
    Returns
    -------
    None
    """
    

    
    
    print('\n' + '*' * (50))
    print('*' + ' ' * (48) + '*')
    print('*' + ' ' * ((48 - 12) // 2) + 'T-one-T-tune' + ' ' * ((48 - 12) // 2) + '*')
    print('*' + ' ' * (48) + '*')
    if module:
        length = len(module) + 7 # length of the module name + " module"
        isdivisible = (48 - length) % 2 == 0
        if isdivisible:
            print('*' + ' ' * ((48 - length) // 2) + f'{module} module' + ' ' * ((48 - length) // 2) + '*')
        else:
            print('*' + ' ' * ((48 - length) // 2) + f'{module} module' + ' ' * ((48 - length) // 2 + 1) + '*')
    print('*' + ' ' * (48) + '*')
    print('**************************************************')
    print('* Authors:                                       *')
    print('*     Enrico Ravera, University of Florence      *')
    print('*     Francesco Bruno, University of Florence    *')
    print('*     Letizia Fiorucci, MPI KoFo                 *')            
    print('**************************************************\n')

def the_end(CO):
    """
    Display the references for the methods used in this software.
    
    Parameters
    ----------
    CO : Conf_Optns object
        The configuration object containing the references and DOIs.
    
    Returns
    -------
    None
    """
    
    print('\n*****************************************************\n')
    print('Analysis complete. Below are the references for the methods used in this software.')
    print('We can not nor we want to force you to cite any of these papers, but we - and likewise (we believe) the authors of the original papers - appreciate if you do.')

    for ref, doi in zip(CO.citelist, CO.doilist):
        print(f'- {ref} DOI: {doi}\n') 

    print('*****************************************************\n')    
    
class Conf_Optns:
    refdictionary = {'fushman': ['Fushman, D. (2012). Determining Protein Dynamics from 15N Relaxation Data by Using DYNAMICS. In: Shekhtman, A., Burz, D. (eds) Protein NMR Techniques. Methods in Molecular Biology, vol 831. Humana Press.','10.1007/978-1-61779-480-3_24'],
                    'cavanagh': ['Cavanagh, J., Fairbrother, W. J., Palmer III, A. G., Rance, M., & Skelton, N. J. (2007). Protein NMR spectroscopy: principles and practice. Academic Press.','10.1016/B978-0-12-164491-8.X5000-3'],
                    'bermel': ['Bermel, W., Bertini, I., Felli, I. C., Piccioli, M., Pierattelli, R. (2006). Progress in NMR Spectroscopy, 48, 25-45','doi.org/10.1016/j.pnmrs.2005.09.002'],
                    'bolognesi': ['Bolognesi, T., Schiavina, M., Felli, I. C., & Pierattelli, R. (2025). NMR insights on multidomain proteins: the case of the SARS-CoV-2 nucleoprotein. Progress in Nuclear Magnetic Resonance Spectroscopy, 148, 101577.','10.1016/j.pnmrs.2025.101577'],
                    'ferrage': ['Ferrage, F. (2011). Protein dynamics by 15N nuclear magnetic relaxation. In: Shekhtman, A., Burz, D. (eds) Protein NMR Techniques. Methods in Molecular Biology, vol 831. Humana Press.','10.1007/978-1-61779-480-3_9'],
                    'lee': ['Lee, D., Hilty, C., Wider, G., Wüthrich, K. (2006). J. Magn. Reson., 178, 72-76','10.1016/j.jmr.2005.08.014'],
                    'robson': ['Robson, S. A., Dağ, Ç., Wu, H., & Ziarek, J. J. (2021). TRACT revisited: an algebraic solution for determining overall rotational correlation times from cross-correlated relaxation rates. Journal of biomolecular NMR, 75(8), 293-302.','10.1007/s10858-021-00379-5'],
                    'Fiorucci': ['Fiorucci, L., Bruno, F., et al. (2024). TrAGICo: a tool for the analysis of TRACT experiments and the design of T1 and T2 experiments. Magnetic Resonance in Chemistry, 62(6), e5537.','10.1002/mrc.5537'],
                    'rezaei-ghaleh' : ['Rezaei‐Ghaleh, N., Parigi, G., Soranno, A., Holla, A., Becker, S., Schuler, B., ... & Zweckstetter, M. (2018). Local and global dynamics in intrinsically disordered synuclein. Angewandte Chemie International Edition, 57(46), 15262-15266.', '10.1002/anie.201808172'],
                    'salvi' : ['Salvi, N., Abyzov, A., Blackledge, M. (2017). Atomic resolution conformational dynamics of intrinsically disordered proteins from NMR spin relaxation. Progress in Nuclear Magnetic Resonance Spectroscopy, 102-103, 43-60.', '10.1016/j.pnmrs.2017.06.001'],
                    'parigi2000' : ['Bertini, I., Fragai, M., Luchinat, C., & Parigi, G. (2000). 1H NMRD profiles of diamagnetic proteins: a model‐free analysis. Magnetic Resonance in Chemistry, 38(7), 543-550.', '10.1002/1097-458X(200007)38:7<543::AID-MRC722>3.0.CO;2-%23'],
                    'klassez': ['Bruno, F. (2025). Automatized quantitative NMR for industrial applications. Doctoral Thesis, University of Florence', 'None'],
                     }
    def __init__(self, parser, module='makelists'):
        """
        Initialize the CO object with the command-line arguments and evaluate the true/false options.
        
        Parameters
        ----------
        parser : argparse.Namespace
            The parser object containing the command-line arguments.
        
        Returns
        -------
        None
            The function initializes the CO object with the command-line arguments
        """

        self.module = module
        splashscreen(module=self.module)
        self.eval_truefalse(parser)
        self.evaluate_S2_tau(parser)
        self.check_values(parser)
        self.get_experiment(parser, config_p=load_config())
        
        
    def add_ref(self, ref):
        """
        Adds a reference to the list of references to be cited. The reference is added to the ``citelist`` attribute of the ``CO`` object and its DOI is added to the ``doilist`` attribute. The reference is identified by its key in the ``refdictionary`` attribute of the ``CO`` object, which contains the full reference and its DOI.

        Parameters
        ----------
        ref : str
            The key of the reference to be added, which should be present in the ``refdictionary`` attribute of the ``CO`` object.
            
        Returns
        -------
        None
            The function updates the ``citelist`` and ``doilist`` attributes of the ``CO`` object with the reference and its DOI, respectively. If the reference is already present in the ``doilist``, it is not added again.
                
        """
        if not hasattr(self, 'doilist'):
            self.doilist = []
        if not hasattr(self, 'citelist'):
            self.citelist = []
        if self.refdictionary[ref][1] not in self.doilist:
            self.doilist.append(self.refdictionary[ref][1])
            self.citelist.append(self.refdictionary[ref][0])

    def eval_truefalse(self, parser):
        """
        Evaluate the true/false options and store them in the ``options`` attribute of the ``CO`` object. The options are:
        
            - integrate: whether to perform the analysis using integrals of the spectra instead of point-by-point fitting. This option is used in the analysis of TRACT experiments.
            - selectregion: whether to select a region of the spectrum for the analysis instead of using a default range. This option is used in the analysis of TRACT experiments.
            - phase: whether to phase the spectra. This option is used in the analysis of TRACT experiments and for the NS module.
            - smoothdata: whether to apply a Savitzky-Golay filter to the spectra before the analysis. This option is used in the analysis of TRACT experiments.
            - smoothrates: whether to apply a Savitzky-Golay filter to the relaxation rates before the design of the T1 and T2 experiments. This option is used in the design of T1 and T2 experiments.
            - plot: whether to plot the spectra and the results of the analysis. This option is used in the analysis of TRACT experiments.
            - logscale: whether to create a logarithmically spaced vdlist for the design of T1 and T2 experiments. This option is used in the design of T1 and T2 experiments and in the setup of TRACT experiments.
            - large: whether to create the lists for the "large" sequence, which is optimized for short T2 times. This option is used in the design of T1 and T2 experiments. If True, the d21 value is set to 450 us and only 8 cycles per CPMG block are used instead of 16.
            - small: whether to create the lists for the "small" sequence, which is optimized for long T2 times. This option is used in the design of T1 and T2 experiments. If True, the d21 value is set to 600 ms and 16 cycles per CPMG block are used.
            - idp: whether to use the IDP model to extract the order parameter S2 instead of tau_c. This option is used in the analysis of TRACT experiments and in the design of T1 and T2 experiments. If True, the molecular weight of the protein must be provided using the --MW argument.
            - randomize: whether to randomize the order of the values in the vdlist for the design of T1 and T2 experiments. This option is used in the design of T1 and T2 experiments.

        Parameters
        ----------
        parser : argparse.Namespace
            The parser object containing the command-line arguments.

        Returns
        -------
        None
            Updates the ``options`` attribute of the ``CO`` object with the evaluated options and adds the references. If the options are mutually exclusive, it raises a ValueError with an appropriate message.
        """
        self.options = {
            'integrate': parser.integrate if hasattr(parser, 'integrate') else None, # tract options
            'selectregion': parser.selectregion if hasattr(parser, 'selectregion') else None, #tract options
            'phase': parser.phase if hasattr(parser, 'phase') else None, #tract and ns options
            'smoothdata': parser.smoothdata if hasattr(parser, 'smoothdata') else None, # tract options
            'smoothrates': parser.smoothrates if hasattr(parser, 'smoothrates') else None, # tract options
            'plot': parser.plot if hasattr(parser, 'plot') else None, #tract options
            'logscale': parser.logscale if hasattr(parser, 'logscale') else None, # makelists or setuptract options
            'large': parser.large if hasattr(parser, 'large') else None, # makelists option
            'small': parser.small if hasattr(parser, 'small') else None, # makelists option
            'idp': parser.idp if hasattr(parser, 'idp') else None, # all except ns
            'randomize': parser.randomize if hasattr(parser, 'randomize') else None, # makelists option
            'readints': parser.readints if hasattr(parser, 'readints') else None # tract options
        }
        if self.module != 'tract':
            if self.options['integrate'] or self.options['selectregion']:
                raise ValueError('The --integrate, --selectregion and --phase options are only valid for the TRACT module. Please remove them and run the script again.')
        else:
            if self.options['smoothdata'] and self.options['smoothrates']:
                raise ValueError('The --smoothdata and --smoothrates options are mutually exclusive, as they are used to apply a Savitzky-Golay filter to the spectra or to the relaxation rates, respectively. Please choose one of the two options and run the script again.')
            if self.options['integrate'] and self.options['selectregion']:
                raise ValueError('The --integrate and --selectregion options are mutually exclusive, as they are used to apply different methods for the analysis of TRACT experiments. Please choose one of the two options and run the script again.')
            if self.options['integrate'] is None and self.options['selectregion'] is None:
                print('Neither --integrate nor --selectregion options were provided. Defaulting to --selectregion.')
                self.options['integrate'] = False
                self.options['selectregion'] = True   
        if self.module != 'makelists':
            if self.options['large'] or self.options['small']:
                raise ValueError('The --large and --small options are only valid for the makelists module. Please remove them and run the script again.')
        else:
            if self.options['large'] and self.options['small']:
                raise ValueError('The --large and --small options are mutually exclusive, as they are used to suggest different vdlist for the design of T1 and T2 experiments. Please choose one of the two options and run the script again.')        
        
        
    def evaluate_S2_tau(self, parser):
        """
        Evaluate the order parameter S2 and the correlation time tau to be used for the calculations.
        If the --idp option is used, check if the molecular weight of the protein is provided using the --MW argument. 
        If not, raise a ValueError with an appropriate message. I
        If the molecular weight is provided, calculate the tau_slow and store it in the first entry of ``tau`` attribute of the ``CO`` object. 
        Also, add the reference for the IDP model to the list of references to be cited. 
        If the --S2 option is provided, in the TRACT modules use this value for the order parameter S2 instead of the default values of 0.15. In the other modules 2 values should be provided with the --idp flag. They are then stored in the ``S2`` attribute of the ``CO`` object. 
        If the --corr_window_idp option is provided, calculate the correlation time of the intermediate motion using the formula MWi = corr_window_idp * 0.110 and computes the intermediate tau from it. Otherwise the tau defaults to 1.6 ns.

        Parameters
        ----------
        parser : argparse.Namespace
            The parser object containing the command-line arguments.
        Returns
        -------
        None
            The function updates the ``S2`` and ``tau`` attributes of the ``CO`` object, and adds the references to the list of references to be cited.
        
        """       
        self.S2 = [0.9]        
        if self.options['idp']:
            if not hasattr(parser, 'MW') and not hasattr(parser, 'tau') or (parser.MW is None and parser.tau is None):
                raise ValueError('Molecular weight or tau of the protein not provided. Please provide it using the --MW or --tau argument.')
            else:
                self.S2 = [0.15, 0.3]
                if parser.S2:
                    self.S2[0] = float(parser.S2[0])
                    if len(parser.S2) > 1:
                        self.S2[1] = float(parser.S2[1])
                    self.add_ref('rezaei-ghaleh')
            if parser.corr_window_idp:
                self.MWi = float(parser.corr_window_idp) * 0.110
        else:
            if parser.S2:
                self.S2[0] = float(parser.S2[0])
                
        self.T = 298.15
        if parser.T:
            self.T = float(parser.T)
        if hasattr(parser, 'tau') and parser.tau is not None:
            self.tau = [float(t)*1e-9 for t in parser.tau]
        if hasattr(parser, 'MW') and parser.MW is not None:
            self.MW = float(parser.MW)
            if not hasattr(self, 'tau'):
                self.tau = [(self.MW * 0.5998 + 0.1674) * 1e-9]
                self.add_ref('cavanagh')
                if self.T != 298.15:
                    self.tau[0] = hydrodynamics_utils.recompute_tau(self.tau[0], self.T)
        if self.options['idp']:
            if not hasattr(self, 'tau'):
                self.tau = [(self.MW * 0.5998 + 0.1674) * 1e-9]
                self.add_ref('cavanagh')
                if self.T != 298.15:
                    self.tau[0] = hydrodynamics_utils.recompute_tau(self.tau[0], self.T)
            if len(self.tau) == 1:
                self.tau.append((self.MWi * 0.5998 + 0.1674) * 1e-9 if hasattr(self, 'MWi') else 1.6e-9) # default value for the correlation time of the intermediate motion in the IDP model is 1.6 ns
                self.add_ref('rezaei-ghaleh')
        if hasattr(self, 'tau'):
            self.tau.append(1e-11) # the fast motion is always set to 10 ps.
        if hasattr(parser, 'basedir') and parser.basedir is not None:
            self.basedir = parser.basedir[0]
        if hasattr(parser, 'tract') and parser.tract is not None:
            self.tract = parser.tract[0]
        
    def check_values(self, parser):
        """
        Check the values provided by the user for compatibility with the options selected.

        Parameters
        ----------
        parser : argparse.Namespace
            The parser object containing the command-line arguments.
        Returns
        -------
        None
            The function checks the values provided by the user for compatibility with the options selected and raises ValueError with an appropriate message if any incompatibility is found.
            Updates the attributes of the ``CO`` object with the values provided by the user.
        """
        
        if hasattr(parser, 'nT'):
            if self.module == 'setuptract':
                if parser.nT is None:
                    self.nT = 8
                else:
                    if len(parser.nT) != 1:
                        raise ValueError('The --nT option for the setuptract module can only accept one value, corresponding to the number of T1 and T2 increments to be designed. Please provide one value to use the same number of increments for both T1 and T2.')
                    self.nT = [int(parser.nT[0])]
            elif self.module == 'makelists':
                if parser.nT is None:
                    self.nT = [8, 8]
                else:
                    if len(parser.nT) ==1:
                        self.nT = [int(parser.nT[0]), int(parser.nT[0])]
                    elif len(parser.nT) ==2:
                        self.nT = [int(parser.nT[0]), int(parser.nT[1])]
                    else:
                        raise ValueError('The --nT option can only accept one or two values, corresponding to the number of T1 and T2 increments to be designed, respectively. Please provide one value to use the same number of increments for both T1 and T2, or two values to specify different numbers for T1 and T2.')
            elif self.module == 'setupsolventpre':
                if parser.nT is None:
                    self.nT = [8, 2]
                else:    
                    if len(parser.nT) == 1:
                        self.nT = [int(parser.nT[0]), 2]
                    elif len(parser.nT) == 2:
                        self.nT = [int(parser.nT[0]), int(parser.nT[1])]
                    else:
                        raise ValueError('The --nT option can only accept one or two values, corresponding to the number of T1 and T2 increments to be designed, respectively. Please provide one value to use the same number of increments for both T1 and T2, or two values to specify different numbers for T1 and T2.')
            else:
                raise ValueError('The --nT option is only valid for the setuptract, makelists, and setupsolventpre modules. Please remove it and run the script again.')
        if self.options['smoothdata']: 
            if not hasattr(parser, 'slw') or not parser.slw:
                self.slw = 5 # default value for the percentage of the spectrum to use for the sliding window smoothing
            else:
                self.slw = parser.slw[0]
        if self.options['smoothrates']: 
            if not hasattr(parser, 'slw') or not parser.slw:
                self.slw = 5 # default value for the percentage of the spectrum to use for the sliding window smoothing
            else:
                self.slw = parser.slw[0]
        if parser.nucs is not None:
            self.nucs = parser.nucs
        else:
            self.nucs = ['1H', '15N']
        if parser.r is not None:
            self.r = float(parser.r)*1e-10
        else:
            if '1H' in self.nucs and '15N' in self.nucs:
                self.r = 1.02e-10 # in meters, default value for the distance between the 1H and 15N nuclei in a protein amide group
            elif '1H' in self.nucs and '13C' in self.nucs:
                self.r = 1.23e-10 # in meters, default value for the distance between the 1H and 13C nuclei in a protein alpha carbon 
            else:
                raise ValueError('Unsupported combination of nuclei. Please provide the r to be used.')
        if parser.Deltasigma is not None:
            self.Deltasigma = float(parser.Deltasigma)
        else:
            if '15N' in self.nucs and '1H' in self.nucs:
                self.Deltasigma = -160 # default value for the chemical shift anisotropy of the 15N nucleus in a protein amide group
            elif '13C' in self.nucs and '1H' in self.nucs:
                self.Deltasigma = 60 # default value for the chemical shift anisotropy of the 13C nucleus in a protein alpha carbon
            else:
                raise ValueError('Unsupported combination of nuclei. Please provide the Deltasigma to be used.')
        if parser.theta is not None:
            self.theta = float(parser.theta) * np.pi / 180 # convert from degrees to radians
        else:
            if '15N' in self.nucs and '1H' in self.nucs:
                self.theta = 17 * np.pi / 180 # in radians, default value for the angle between the N-H bond and the principal axis of the chemical shift tensor of the 15N nucleus in a protein amide group 17°
            elif '13C' in self.nucs and '1H' in self.nucs:
                self.theta = 109 * np.pi / 180 # in radians, default value for the angle between the C-H bond and the principal axis of the chemical shift tensor of the 13C nucleus in a protein alpha carbon 109°
            else:
                raise ValueError('Unsupported combination of nuclei. Please provide the theta to be used.')
        if hasattr(parser, 'xred') and parser.xred is not None:
            self.xred = [float(x) for x in parser.xred]
        if self.module=='NS':
            if not hasattr(self, 'xred'):
                print('No value provided for the xred parameter. Defaulting to [0.1, 0.3, 0.7].')
                self.xred = [0.1, 0.3, 0.7]
            if parser.lw is not None:
                self.lw = parser.lw
            if parser.nres is not None:
                self.nres = parser.nres
        if self.module!='tract':
            if hasattr(parser, 'B0') and parser.B0 is not None:
                if len(parser.B0) == 1:
                    self.B_0 = float(parser.B0[0])
            
                else:
                    raise NotImplementedError('Multiple values for B0 are not supported yet. Please provide one value for the magnetic field strength to be used in the calculations.')
            else:
                self.get_B0(config_p=load_config())
        
    def get_B0(self, config_p=None):
        """
        Get the magnetic field value from the config file if it exists, otherwise ask the user to provide it.
        The magnetic field value is stored in the ``B_0`` attribute of the ``CO`` object.
        
        Parameters
        ----------
        config_p : configparser.ConfigParser, optional
            The config object containing the parameters. If None, the function will ask the user to provide the magnetic field value. Default is None.
            
        Returns
        -------
        None
            The function updates the ``B_0`` attribute of the ``CO`` object with the magnetic field value in Tesla, either from the config file or from the user input.
        """
        if config_p is not None:
            self.B_0 = float(config_p['DEFAULT']['B0'])        
        if not hasattr(self, 'B_0'):
            self.B_0 = float(input('Please provide the magnetic field value in Tesla: '))
    
    def get_experiment(self, parser, config_p=None):
        """
        Configures the experiment folders for the default experiments.
        
        Parameters
        ----------
        config_p : configparser.ConfigParser, optional
            The config object containing the parameters. If None, the function will ask the user to provide the experiment parameters. Default is None.
            
        Returns
        -------
        None
            The function updates the ``basedir``, ``tract``, and ``hsqc`` attributes of the ``CO`` object with the experiment parameters, either from the config file or from the user input.
        """

        if config_p is not None:
            if 'EXPERIMENT' in config_p.keys():
                self.basedir = config_p['EXPERIMENT']['basedir']
                self.tract = config_p['EXPERIMENT']['expno']
                self.hsqc = config_p['EXPERIMENT']['hsqcexpno']
                print(f'Experiment parameters loaded from config file: basedir={self.basedir}, tract={self.tract}, hsqc={self.hsqc}')
            else:
                print('No EXPERIMENT section found in config file. Please provide the experiment parameters manually.')
                self.basedir = input('Please provide the base directory of the experiment: ')
                self.tract = input('Please provide the experiment number of the TRACT experiment: ')
                self.hsqc = input('Please provide the experiment number of the reference HSQC spectrum: ')
        elif os.getlogin() == 'ravera' and socket.gethostname() == 'dual':
            if hasattr(self, 'MW') and self.MW is not None:
                self.basedir = None
                self.tract = None
                self.hsqc = None
            elif hasattr(self, 'tau') and self.tau is not None:
                self.basedir = None
                self.tract = None
                self.hsqc = None
            else:
                self.basedir = 'tract'
                self.tract = '36'
                self.hsqc = '411'
                print(f'Experiment parameters set for internal test environment: basedir={self.basedir}, tract={self.tract}, hsqc={self.hsqc}')
        else:
            if hasattr(self, 'MW') and self.MW is not None:
                self.basedir = None
                self.tract = None
                self.hsqc = None
            elif hasattr(self, 'tau') and self.tau is not None:
                self.basedir = None
                self.tract = None
                self.hsqc = None
            else:
                pass
        if hasattr(parser, 'basedir') and parser.basedir is not None:
            self.basedir = parser.basedir
        if hasattr(parser, 'tract') and parser.tract is not None:
            self.tract = int(parser.tract)
        if hasattr(parser, 'hsqc') and parser.hsqc is not None:
            self.hsqc = int(parser.hsqc)
        