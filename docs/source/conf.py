#! /usr/bin/env python3

import os
import sys
from pathlib import Path

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'T1Ttune'
copyright = '2026, Francesco Bruno, Letizia Fiorucci, Enrico Ravera'
author = 'Francesco Bruno, Letizia Fiorucci, Enrico Ravera'
release = '0.0.1'

package_path = Path(__file__).parent.parent.parent.absolute()
print('!!!', package_path)
sys.path.insert(0, package_path)  # Adjust the path to include the parent directory

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = []

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
}

autosummary_generate = True
autoclass_content = 'class'
autodoc_member_order = 'alphabetical'
autodoc_default_options = {
        'members': True,
        'undoc_members': True,
        'show_inheritance': True,
        'special_members': '__init__,__call__',
        'member_order': 'alphabetical',
        }


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for LaTeX output ------------------------------------------------

# Options for LaTeX output
# ------------------------

# The paper size ('letter' or 'a4').
latex_paper_size = 'letter'

# Grouping the document tree into LaTeX files.
# List of tuples:
#   (source start file, target name, title, author,
#    document class [howto/manual])

latex_documents = [
    ('index', 't1ttune.tex', 'T1Ttune Documentation',
     'Enrico Ravera, Francesco Bruno, and Letizia Fiorucci', 'manual'),
]


# The name of an image file (relative to this directory) to place at the top of
# the title page.
latex_logo = None

# Use Unicode aware LaTeX engine
latex_engine = 'pdflatex'  # or 'lualatex'

latex_elements = {}

# Keep babel usage also with xelatex (Sphinx default is polyglossia)
# If this key is removed or changed, latex build directory must be cleaned
latex_elements['babel'] = r'\usepackage{babel}'

# Font configuration
# Fix fontspec converting " into right curly quotes in PDF
# cf https://github.com/sphinx-doc/sphinx/pull/6888/

# Fix fancyhdr complaining about \headheight being too small
latex_elements['passoptionstopackages'] = r"""
    \PassOptionsToPackage{headheight=14pt}{geometry}
"""

# Additional stuff for the LaTeX preamble.
latex_elements['preamble'] = r"""
   % Show Parts and Chapters in Table of Contents
   \setcounter{tocdepth}{0}
   % One line per author on title page
   \DeclareRobustCommand{\and}%
     {\end{tabular}\kern-\tabcolsep\\\begin{tabular}[t]{c}}%
   \usepackage{etoolbox}
   \let\latexdescription=\description
   \def\description{\latexdescription{}{} \breaklabel}
   % But expdlist old LaTeX package requires fixes:
   % 1) remove extra space
   \makeatletter
   \patchcmd\@item{{\@breaklabel} }{{\@breaklabel}}{}{}
   \makeatother
   % 2) fix bug in expdlist's way of breaking the line after long item label
   \makeatletter
   \def\breaklabel{%
       \def\@breaklabel{%
           \leavevmode\par
           % now a hack because Sphinx inserts \leavevmode after term node
           \def\leavevmode{\def\leavevmode{\unhbox\voidb@x}}%
      }%
   }
   \newcommand{\ui}{\mathrm{i}}
   \makeatother
"""
# Sphinx 1.5 provides this to avoid "too deeply nested" LaTeX error
# and usage of "enumitem" LaTeX package is unneeded.
# Value can be increased but do not set it to something such as 2048
# which needlessly would trigger creation of thousands of TeX macros
latex_elements['maxlistdepth'] = '10'
latex_elements['pointsize'] = '11pt'

# Better looking general index in PDF
latex_elements['printindex'] = r'\footnotesize\raggedright\printindex'

# Documents to append as an appendix to all manuals.
latex_appendices = []

# If false, no module index is generated.
latex_use_modindex = True

latex_toplevel_sectioning = 'chapter'



autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,
    'show-inheritance': True,
}

numpydoc_show_class_members = False
numpydoc_xref_param_type = True

