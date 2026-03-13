#! /usr/bin/env python3

"""
Registry dei sottocomandi.

"""

from t1ttune.scripts.t1ttune_configure import ConfigureCmd
from t1ttune.scripts.t1ttune_makelists import MakeListsCmd
from t1ttune.scripts.t1ttune_setuptract import SetupTractCmd
from t1ttune.scripts.t1ttune_tract import TractCmd
from t1ttune.scripts.t1ttune_ns import NSCmd
from t1ttune.scripts.t1ttune_interactive import InteractiveCmd
from t1ttune.scripts.t1ttune_solventpre import SPREListsCmd

registry = {
    "makelists": MakeListsCmd,
    "solventpre": SPREListsCmd,
    #"shuttle": ShuttleCmd,
    "setuptract": SetupTractCmd,    
    "configure": ConfigureCmd,
    "ns": NSCmd,    
    "tract": TractCmd,
    "interactive": InteractiveCmd,
}

