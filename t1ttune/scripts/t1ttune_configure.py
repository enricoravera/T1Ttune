#! /usr/bin/env python3

import os
import configparser
import klassez as kz
import argparse

from .base import BaseCommand

from .t1ttune_utils import istract
from . import f_findfs


class ConfigureCmd(BaseCommand):
    SHORT_HELP = "Configure T1tune for the current system"
    DESCRIPTION = '\n'.join([
        "This command checks if the software is installed on a spectrometer workstation.",
        "If it is, it creates a configuration file with the appropriate paths and settings for that system."
    ])
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--topspinpath', type=str, help='the path to the topspin installation. If not provided, the script will try to find it automatically and use the found path as default value for this argument.')
        parser.add_argument('--basedir', type=str, help='the directory of the default example experiments.')
        parser.add_argument('--hsqc', type=str, help='if set, the command will also ask for the expno of the HSQC experiment in the standard dataset.')
        parser.add_argument('--tract', type=str, help='if set, the command will also ask for the expno of the TRACT experiment in the standard dataset.')
    @staticmethod
    def run(args: argparse.Namespace) -> None:
        if args.topspinpath:
            fspath = args.topspinpath
            if not os.path.isdir(fspath):
                raise ValueError(f"The provided path {fspath} is not a valid directory.")
            if not os.path.isdir(os.path.join(fspath, "conf", "instr", "spect")):
                if not os.path.exists(os.path.join(fspath, "conf", "instr", "remote_spect", "uxnmr.info")):
                    raise ValueError(f"The provided path {fspath} does not appear to be a valid TopSpin installation on a spectrometer workstation.")
        else:
            findfs_info = f_findfs.find_topspin()
            if not findfs_info["found"]:
                raise RuntimeError("TopSpin does not appear to be installed on this workstation. Cannot configure T1tune.")
            if not findfs_info["spectrometer"]:
                raise RuntimeError("TopSpin does not appear to be installed on a spectrometer workstation. Cannot configure T1tune.")
            fspath = findfs_info["install_path"]
            
        config(fspath, args.basedir, args.hsqc, args.tract)
        
def config(foundpath, basedir=None, hsqcexpno=None, tractexpno=None):
    """
    Configures the script for use on a spectrometer, and is reserved to the NMR superuser.
    The function prompts the user to input the path to the topspin installation and checks if the user is in the list of NMR superusers. 
    If the user is not a superuser, a PermissionError is raised and the function exits. 
    If the user is a superuser, the function prompts the user to input the path to a standard dataset, the experiment number of the TRACT experiment in the standard dataset, and the experiment number of the HSQC experiment in the standard dataset (used only to suggest a number of scans for the experiments). 
    The function then loads the TRACT experiment from the standard dataset, checks if it is a TRACT experiment, and extracts B0 field strength from the acqus file. 
    Finally, it creates a config file in the curdir folder of topspin with the extracted information and exits.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
        The function creates a config file in the curdir folder of topspin with the extracted information and exits.
    """
    
    fspath = input('Please provide the path to topspin: ') or foundpath 
    nmrsupath = os.path.join(fspath, 'conf', 'nmrsuperuser')
    with open(nmrsupath) as file:
        nmrsulist = [line.rstrip() for line in file]
    #check if the user corresponds to the admin
    if os.getlogin() not in nmrsulist:
        raise PermissionError('Only the NMR superuser can write the config file. Exiting.')
    else:
        config_p = configparser.ConfigParser()
        if not os.path.isdir(os.path.join(fspath, "conf", "instr", "spect")):
            if not os.path.exists(os.path.join(fspath, "conf", "instr", "remote_spect", "uxnmr.info")):
                raise NameError(f"The provided path {fspath} does not appear to be a valid TopSpin installation on a spectrometer workstation.")
            else:
                spect_folder = os.path.join(fspath, "conf", "instr", "remote_spect") # remote instrument config
        else:
            spect_folder = os.path.join(fspath, "conf", "instr", "spect") # instrument config
        with open(os.path.join(spect_folder, "uxnmr.info")) as spectfile:
            spectinfolist = spectfile.readlines()
            for line in spectinfolist:
                if "1H-frequency" in line:
                    larmor_freq = float(line.split(":")[1].strip(" MHz\n"))
                    break
        B_0_config = larmor_freq / kz.sim.gamma["1H"] # in Tesla, calculated from the 1H Larmor frequency
        config_p['DEFAULT'] = {
            'B0': B_0_config            
        }
        pathtodataset = input('Please provide the path to a standard dataset to be used as example') or basedir
        if not os.path.isdir(pathtodataset):
            raise ValueError(f"The provided path {pathtodataset} is not a valid directory.")
        tractexpno = input('Please provide the experiment number of the TRACT experiment in the standard dataset') or tractexpno
        hsqcexpno = input('Please provide the experiment number of the HSQC experiment in the standard dataset (used only to suggest a number of scans for the experiments)') or hsqcexpno
        if tractexpno is None and hsqcexpno is None:
            raise ValueError('At least one of the two arguments --tract and --hsqc must be provided to specify the experiment numbers of the TRACT and HSQC experiments in the standard dataset, if the standard dataset is provided.')

        S_config= kz.Pseudo_2D(os.path.join(pathtodataset, f'{tractexpno}'))
        if not istract(S_config):
            raise NameError(f'Experiment {tractexpno} is not a TRACT experiment')

        config_p['EXPERIMENT'] = {
            'basedir': pathtodataset,
            'expno': tractexpno,
            'hsqcexpno': hsqcexpno,
        }
        pathtocurdir = os.path.join(fspath, 'prog', 'curdir', os.getlogin(), 'tract_analysis_config.ini')
        with open(pathtocurdir, 'w') as configfile:
            config_p.write(configfile)
    exit()

