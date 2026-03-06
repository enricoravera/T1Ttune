#! /usr/bin/env python3

import os
import sys
import platform
from pathlib import Path
import glob
from itertools import chain
from packaging import version

from .textcolor import textcolor


def find_topspin() -> dict:
    """
    Detect whether Bruker TopSpin is installed on the current workstation and if it is a spectrometer installation.
    Works cross-platform: Windows, Linux, macOS.

    Returns a dict with:
        - found (bool): whether TopSpin was detected
        - version (str|None): version string if detectable
        - install_path (Path|None): path to the TopSpin installation
        - spectrometer (bool): whether this installation appears to be on a spectrometer workstation (based on presence of conf/instr/spect)
    """
    result = {"found": False, "version": None, "install_path": None, "spectrometer": False}

    path = None
    system = platform.system()

    # ------------------------------------------------------------------
    # 1. Check well-known default installation paths
    # ------------------------------------------------------------------
    candidate_paths: list[Path] = []

    if system == "Windows":
        # Default Bruker install roots on Windows
        for drive in ["C:", "D:"]:
            candidate_paths += [
                Path(drive) / "Bruker" / "TopSpin",
                Path(drive) / "TopSpin",
            ]
        # Also check 32-bit Program Files on 64-bit Windows
        for pf in ["ProgramFiles", "ProgramFiles(x86)", "ProgramW6432"]:
            pf_dir = os.environ.get(pf)
            if pf_dir:
                candidate_paths += [
                    Path(pf_dir) / "Bruker" / "TopSpin",
                ]

    elif system == "Linux":
        candidate_paths += [
            Path("/opt/topspin"),
            Path("/opt/Bruker/topspin"),
            Path.home() / "topspin",
        ]

    elif system == "Darwin": # macOS (rare but possible)
        candidate_paths += [
            Path("/Applications/Bruker/TopSpin"),
            Path("/Applications/TopSpin"),
            Path.home() / "Applications" / "TopSpin",
        ]


    inst_folder = os.path.join("exp", "stan", "nmr")
    spect_folder = os.path.join("conf", "instr", "spect") # instrument config
    uxnmrinfofile = os.path.join("conf", "instr", "remote_spect", "uxnmr.info") # remote instrument config

    globbed = []
    for p in candidate_paths:
        globbed.append(glob.glob(str(p) + '*'))
    globbed = list(set(list(chain.from_iterable(globbed))))

    if len(globbed) == 0:
        return result
    
    versions = []
                
    for i, g in enumerate(globbed):
        if (Path(g) / inst_folder).exists():
            v = version.parse(g.rsplit("topspin", 1)[-1])
            if (Path(g) / spect_folder).exists():
                is_spectrometer = True
            elif (Path(g) / uxnmrinfofile ).exists():
                is_spectrometer = True
            else:
                is_spectrometer = False
            versions.append([v, g, is_spectrometer])
        else:
            continue
    if len(versions) == 0:
        return result
    else:
        versions.sort(key=lambda x: x[0], reverse=True) # sort by version, descending
        for v, g, is_spectrometer in versions:
            if is_spectrometer:
                path = g
                latest_version = v
                break
        else:            # If no spectrometer installation found, take the latest non-spectrometer one
            latest_version, path, is_spectrometer = versions[0]
    if path is not None:
        result.update(
            found=True,
            install_path=path,
            version=str(latest_version),
            spectrometer=is_spectrometer,
        )
    return result

def version_dict(found=False, version=None, install_path=None, spectrometer=False) -> dict:
    """
    Helper function to create a version dict with the given values, using defaults for any missing values.
    """
    return {
        "found": found,
        "version": version,
        "install_path": install_path,
        "spectrometer": spectrometer,
    }

# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    info = find_topspin()
    if info["found"]:
        print(textcolor("TopSpin is installed.", "green"))
        print(f" Path : {info['install_path']}")
        print(f" Version : {info['version'] or 'unknown'}")
        print(f" Is Spectrometer: {info['spectrometer']}")
    else:
        print(textcolor("TopSpin does not appear to be installed on this workstation.", "orange"))
    sys.exit(0 if info["found"] else 1)
