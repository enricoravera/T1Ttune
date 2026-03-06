#! /usr/bin/env python3

def textcolor(text, color, bold=False):
    """
    Change the color and formatting of the text for terminal output.

    Parameters
    ----------
    text : str
        The text to be colored
    color : str
        The color to apply. Options are ``'red'``, ``'green'``, ``'blue'``, ``'yellow'``, ``'orange'``, ``'magenta'``, ``'default'``.
    bold : bool, optional
        Whether to apply bold formatting. Default is ``False``.

    Returns
    -------
    str
        The colored text.
    """    
    colors = {'red': '\033[38;5;196m',
              'green': '\033[38;5;2m',
              'blue': '\033[38;5;12m',
              'yellow': '\033[38;5;226m',
              'orange': '\033[38;5;208m',
              'magenta': '\033[38;5;13m',
              'default': ''}
    boldcode = '\033[1m' if bold else ''
    resetcode = '\033[0m'
    return f'{boldcode}{colors[color]}{text}{resetcode}'

