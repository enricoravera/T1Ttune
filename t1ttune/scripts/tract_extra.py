#! /usr/bin/env python3

import numpy as np
from copy import deepcopy

def split_tract(S):
    """Split the interleaved TROSY and ANTITROSY datasets of a TRACT experiment.
    
    Parameters
    ----------
    S : Spectrum object
        The Spectrum object containing the interleaved TROSY and ANTITROSY datasets.

    Returns
    -------
    Sa, Sb : Spectrum objects
        The two Spectrum objects containing the TROSY and ANTITROSY datasets, respectively.
    """
    S.fid = np.reshape(S.fid.flatten(), (2*S.fid.shape[0], -1))
    Sa = deepcopy(S)
    Sb = deepcopy(S)
    Sa.fid = S.fid[::2]
    Sb.fid = S.fid[1::2]
    Sa.acqus['TD1'] = Sa.fid.shape[0]
    Sb.acqus['TD1'] = Sb.fid.shape[0]
    return Sa, Sb

