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
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OOFEM Material Library Manual'
copyright = '2026, Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
author = 'Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
release = '3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# doc/_ext holds the extensions shared by the OOFEM manuals; oofemroles defines
# the roles that mirror the manuals' LaTeX macros (\param, \elemparam,
# \descitem, ...).
sys.path.insert(0, os.path.abspath(os.path.join('..', '_ext')))

extensions = ['oofemroles', 'oofemmath', 'oofemlatex']

# The source directory also holds the LaTeX sources and the latex2html output,
# so keep the builder away from everything that is not part of this manual.
exclude_patterns = [
    '_build',
    'html',
    'auto',
    'figs',
    'Thumbs.db',
    '.DS_Store',
    '*.rst.orig',
]

# ":numref:" is used throughout to reference tables and figures by number.
numfig = True

# Math macros of this manual, inherited from matlibmanual.tex.  The shared ones
# (\mbf, \del, \der, \pard) come from oofemmath; declared here are the
# ones specific to this manual.  oofemmath feeds them to MathJax for the
# HTML and to the LaTeX preamble for the PDF, so the two cannot drift.
oofem_math_macros = {
    'alphaPsi':  r'\alpha_{\psi}',
    'bsig':      r'\mbf{\sigma}',
    'dO':        r'\,\mbox{d}\Omega',
    'dvepst':    r'\delta\tilde{\veps}',
    'dvet':      r'\delta\vet',
    'dvs':       r'\delta\vs',
    'dvsig':     r'\delta\vsig',
    'e':         r'\mbf{\varepsilon}',
    'ep':        r'\mbf{\varepsilon}^p',
    'epd':       r'\dot{\mbf{\varepsilon}}^p',
    'eps':       r'\mbf{\varepsilon}',
    'epsp':      r'\eps_{\mathrm{p}}',
    'epspd':     r'\dot{\eps}_{\mathrm{p}}',
    'epss':      r'\varepsilon',
    'fc':        r'\bar{f}_c',
    'ft':        r'\bar{f}_t',
    'kap':       r'\mbf{\kappa}',
    'kappac':    r'\kappa_{\mathrm{c}}',
    'mD':        r'\mbf{D}',
    'qh':        r'q_{\rm h}',
    'quarter':   r'\frac{1}{4}',
    'sig':       r'\mbf{\sigma}',
    'sigs':      r'\sigma',
    'sym':       r'_{\mbox{\small sym}}',
    'tauY':      r'\tau_{\mathrm {Y}}',
    'tenss':     [r'\boldsymbol{#1}', 1],
    'ud':        r'\mathrm{d}',
    've':        r'\mbf{e}',
    'veps':      r'\mbf{\varepsilon}',
    'vepst':     r'\tilde{\veps}',
    'vet':       r'\tilde{\ve}',
    'vs':        r'\mbf{s}',
    'vsig':      r'\mbf{\sigma}',
    'vsigt':     r'\tilde{\vsig}',
    'vst':       r'\mbf{s}^T',
    'vx':        r'\mbf{x}',
    'vxi':       r'\mbf{\xi}',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme, the shared stylesheet and the page width come from global_conf.py
# and doc/_static/oofem.css, so nothing is set here.

master_doc = 'index'
