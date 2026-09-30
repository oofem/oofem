# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

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

extensions = ['oofemroles', 'oofemtikz']

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

# The manual's equations use the macros defined in elementlibmanual.tex,
# so MathJax has to be taught the same definitions.  They are not all
# in the preamble: several are declared part-way through the text.
mathjax3_config = {
    'tex': {
        'macros': {
            'del':       [r'\displaystyle\frac{#1}{#2}', 2],
            'der':       [r'\frac{{\rm d}{#1}}{{\rm d}{#2}}', 2],
            'mbf':       [r'\boldsymbol{#1}', 1],
            'pard':      [r'\frac{\partial{#1}}{\partial{#2}}', 2],
        },
    },
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'
# doc/_static is shared by the OOFEM manuals.
html_static_path = ['../_static']
html_css_files = ['oofem.css']

# The element summary tables are wide, so use the whole browser window.
# alabaster caps the outer container at page_width (940px by default), which
# overrides body_max_width on its own, so both have to be set.  body_max_width
# stays at 100% so the text fills the page width rather than 80% of it.
html_theme_options = {
    'page_width': '80%',
    'body_max_width': '100%',
    'sidebar_width': '250px',
}

master_doc = 'index'
