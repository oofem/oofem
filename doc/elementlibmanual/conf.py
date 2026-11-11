# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# Import shared settings (theme, stylesheet, numfig); the project
# information below overrides the generic one.
from global_conf import *

# -- Project information -----------------------------------------------------

project = 'OOFEM Element Library Manual'
copyright = '2026, Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
author = 'Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
release = '3.0'

# -- General configuration ---------------------------------------------------

# doc/_ext holds the extensions shared by the OOFEM manuals: oofemroles defines
# the roles mirroring the \param, \field and \descitem macros, and oofemtikz
# provides the "tikz" directive used by the element figures.
sys.path.insert(0, os.path.abspath(os.path.join('..', '_ext')))

extensions = ['oofemroles', 'oofemtikz', 'oofemmath', 'oofemlatex']

# The source directory also holds the LaTeX sources and the latex2html output.
exclude_patterns = [
    '_build',
    'html',
    'auto',
    'Thumbs.db',
    '.DS_Store',
]

# ":numref:" is used throughout to reference tables and figures by number.
numfig = True

# The equations of this manual use \mbf, \del, \der and \pard, inherited from
# elementlibmanual.tex.  All four are common to the OOFEM manuals and come from
# oofemmath, which defines them for MathJax (HTML) and in the LaTeX preamble
# (PDF) alike, so no macros need to be declared here.  Anything specific to this
# manual would go in an "oofem_math_macros" dict.

# -- Options for HTML output -------------------------------------------------

# The theme, the shared stylesheet and the page width come from global_conf.py
# and doc/_static/oofem.css, so nothing is set here.

master_doc = 'index'
